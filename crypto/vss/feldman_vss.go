// Copyright © 2019 Binance
//
// This file is part of Binance. The full Binance copyright notice, including
// terms governing use, modification, and redistribution, is contained in the
// file LICENSE at the root of the source code distribution tree.

// Feldman VSS, based on Paul Feldman, 1987., A practical scheme for non-interactive verifiable secret sharing.
// In Foundations of Computer Science, 1987., 28th Annual Symposium on. IEEE, 427–43
//

package vss

import (
	"errors"
	"fmt"
	"math/big"

	"github.com/binance-chain/tss-lib/common"
	matrix "github.com/binance-chain/tss-lib/common/matrix"
	"github.com/binance-chain/tss-lib/crypto"
	"github.com/binance-chain/tss-lib/tss"
)

type (
	Share struct {
		Threshold int
		ID,       // xi
		Share *big.Int // Sigma i
		Rank uint
	}

	Vs []*crypto.ECPoint // v0..vt

	Shares []*Share
)

// polynomial represents a polynomial of arbitrary degree
type Polynomial struct {
	coefficients []*big.Int
}

type Tag struct {
	XCoord   *big.Int
	DiffTime int
}





var (
	ErrNumSharesBelowThreshold             = fmt.Errorf("not enough shares to satisfy the threshold")
	ErrDIFFERENTIALSHOULDBEPOSITIVEINTEGER = fmt.Errorf("The times of Differential should be a non-negative integer!")
	ErrCANTGETBIRKHOFFCOEFFICIENT          = fmt.Errorf("Can't get BirkhoffCiefficient")

	zero = big.NewInt(0)
	one  = big.NewInt(1)
)

// Returns a new array of secret shares created by Shamir's Secret Sharing Algorithm,
// requiring a minimum number of shares to recreate, of length shares, from the input secret
//
func Create(threshold int, secret *big.Int, indexes []*big.Int, ranks []uint) (Vs, Shares, error) {
	if secret == nil || indexes == nil {
		return nil, nil, fmt.Errorf("vss secret or indexes == nil: %v %v", secret, indexes)
	}
	if threshold < 1 {
		return nil, nil, errors.New("vss threshold < 1")
	}
	num := len(indexes)
	if num < threshold {
		return nil, nil, ErrNumSharesBelowThreshold
	}

	var poly Polynomial
	//poly.coefficients = samplePolynomial(threshold, secret)
	testPoly := make([]*big.Int, 3)
	testPoly[0] = secret
	testPoly[1] = big.NewInt(1)
	testPoly[2] = big.NewInt(1)

	poly.coefficients = testPoly

	poly.coefficients[0] = secret // becomes sigma*G in v
	v := make(Vs, len(poly.coefficients))
	for i, ai := range poly.coefficients {
		v[i] = crypto.ScalarBaseMult(tss.EC(), ai)
	}

	shares := make(Shares, num)
	for i := 0; i < num; i++ {
		/*
			if indexes[i].Cmp(zero) == 0 {
				return nil, nil, fmt.Errorf("party index should not be 0")
			}
		*/

		diffPoly, err := Differential(poly, int(ranks[i]))

		if err != nil {
			return nil, nil, ErrNumSharesBelowThreshold
		}

		share := evaluatePolynomial(threshold, diffPoly.coefficients, indexes[i])
		shares[i] = &Share{Threshold: threshold, ID: indexes[i], Rank: ranks[i], Share: share}
	}
	return v, shares, nil
}

func (share *Share) Verify(threshold int, vs Vs) bool {
	if share.Threshold != threshold || vs == nil {
		return false
	}
	var err error

	var shareTag Tag
	shareTag.DiffTime = int(share.Rank)
	shareTag.XCoord = share.ID

	// Very critical
	degreepoly := len(vs) - 1
	scalarSlice, err := GetLinearEquationCoefficient(shareTag, degreepoly)

	if err != nil {
		return false
	}

	v, err := crypto.ComputeLinearCombinationPoint(scalarSlice, vs, tss.EC())

	if err != nil {
		return false
	}

	/*
		modQ := common.ModInt(tss.EC().Params().N)
		v, t := vs[0], one // YRO : we need to have our accumulator outside of the loop
		for j := 1; j <= threshold; j++ {
			// t = k_i^j
			t = modQ.Mul(t, share.ID)
			// v = v * v_j^t
			vjt := vs[j].SetCurve(tss.EC()).ScalarMult(t)
			v, err = v.SetCurve(tss.EC()).Add(vjt)
			if err != nil {
				return false
			}
		}
	*/

	sigmaGi := crypto.ScalarBaseMult(tss.EC(), share.Share)
	return sigmaGi.Equals(v)
}

func (shares Shares) ReConstruct() (secret *big.Int, err error) {
	
	if shares != nil && shares[0].Threshold > len(shares) {
		return nil, ErrNumSharesBelowThreshold
	}
	modN := common.ModInt(tss.EC().Params().N)

	// x coords
	xs := make([]*big.Int, 0)
	for _, share := range shares {
		xs = append(xs, share.ID)
	}

	secret = big.NewInt(0)
	for i, share := range shares {
		times := one
		for j := 0; j < len(xs); j++ {
			if j == i {
				continue
			}
			sub := modN.Sub(xs[j], share.ID)
			subInv := modN.ModInverse(sub)
			div := modN.Mul(xs[j], subInv)
			times = modN.Mul(times, div)
		}

		fTimes := modN.Mul(share.Share, times)
		secret = modN.Add(secret, fTimes)
	}

	return secret, nil
}

/*

if shares != nil && shares[0].Threshold > len(shares) {
		return nil, ErrNumSharesBelowThreshold
	}

	
	tagList := make([]Tag, len(shares))

	for i := 0; i < len(tagList); i++ {
		tagList[i].XCoord = shares[i].ID
		tagList[i].DiffTime = int(shares[i].Rank)
	}

	
		// x coords
		xs := make([]*big.Int, 0)
		for _, share := range shares {
			xs = append(xs, share.ID)
		}

		secret = zero
		for i, share := range shares {
			times := one
			for j := 0; j < len(xs); j++ {
				if j == i {
					continue
				}
				sub := modN.Sub(xs[j], share.ID)
				subInv := modN.ModInverse(sub)
				div := modN.Mul(xs[j], subInv)
				times = modN.Mul(times, div)
			}

			fTimes := modN.Mul(share.Share, times)
			secret = modN.Add(secret, fTimes)
		}
	

	/*
	birkhoffCoefficient, err := GetBirkhoffCiefficient(tagList, len(shares))

	if err != nil {
		return nil, err
	}

	secret = new(big.Int).Mul(shares[0].Share, birkhoffCoefficient[0])

	for i := 1; i < len(shares); i++ {
		tempValue := new(big.Int).Mul(shares[i].Share, birkhoffCoefficient[i])
		secret.Add(secret, tempValue)
	}
	secret.Mod(secret, tss.EC().Params().N)

	if err != nil {
		return nil, err
	}
	return secret, nil



*/






func samplePolynomial(threshold int, secret *big.Int) []*big.Int {
	q := tss.EC().Params().N
	v := make([]*big.Int, threshold+1)
	v[0] = secret
	for i := 1; i <= threshold; i++ {
		ai := common.GetRandomPositiveInt(q)
		v[i] = ai
	}
	return v
}

// Evauluates a polynomial with coefficients such that:
// evaluatePolynomial([a, b, c, d], x):
// 		returns a + bx + cx^2 + dx^3
//
func evaluatePolynomial(threshold int, v []*big.Int, id *big.Int) (result *big.Int) {
	q := tss.EC().Params().N
	modQ := common.ModInt(q)
	result = new(big.Int).Set(v[0])
	X := big.NewInt(int64(1))
	for i := 1; i < len(v); i++ {
		ai := v[i]
		X = modQ.Mul(X, id)
		aiXi := new(big.Int).Mul(ai, X)
		result = modQ.Add(result, aiXi)
	}
	return
}

///
///
///
///

// GetCoefficientOfLagrange returns the first row of birkhoffMatrix
func GetBirkhoffCiefficient(tagSlice []Tag, threshold int) ([]*big.Int, error) {
	var empty []*big.Int

	// Check these shares has rank enough to recover the secret.
	birkoffInverseMatrix, err := GetBirkhoffMatrix(tagSlice, threshold)

	if err != nil {
		return empty, ErrCANTGETBIRKHOFFCOEFFICIENT
	}
	zeroRow, err := birkoffInverseMatrix.GetRow(0)

	if err != nil {
		return empty, ErrCANTGETBIRKHOFFCOEFFICIENT
	}

	return zeroRow, nil
}

// Generate Birkhoff Matrix accoeding to a tag of shareholders.
func GetBirkhoffMatrix(tagSlice []Tag, threshold int) (*matrix.Matrix, error) {

	matrixBirkhoff, err := generateLinearEquationCoefficientMatrix(tagSlice, threshold)

	if err != nil {
		return nil, err
	}

	birkhoffMatrix, _ := matrix.NewMatrix(tss.EC().Params().N, matrixBirkhoff)
	result, err := birkhoffMatrix.Pseudoinverse()
	fmt.Println(result.GetMatrix()[0][0])

	if err != nil {
		return nil, err
	}

	return result, nil
}

// Establish the coefficient of linear system of Birkhoff systems. The relation of Birkhoff matrix and LinearEquationCoefficientMatrix is
// LinearEquationCoefficientMatrix = the inverse of Birkhoff matrix.
// Assume: share1: DiffTime =0, XCoord=1 share2: XCoord =2, DiffTime = 1, share3: XCoord =3, differTime=2
// Then output is:
// 1^(DiffTime)|_{XCoord}  x^(DiffTime)|_{XCoord} (x^2)^(DiffTime)|_{XCoord}
// [       1                       1                          1   ] DiffTime = 0, XCoord =1
// [       0                       1                          4   ] DiffTime = 1, XCoord =2
// [       0                       0                          2   ] DiffTime = 2, XCoord =3
func generateLinearEquationCoefficientMatrix(tagSlice []Tag, nThreshold int) ([][]*big.Int, error) {

	var result = make([][]*big.Int, len(tagSlice))
	var emptyMatrix = make([][]*big.Int, len(tagSlice))

	if len(tagSlice) < nThreshold {
		return emptyMatrix, ErrCANTGETBIRKHOFFCOEFFICIENT
	}

	for i := int64(0); i < int64(len(tagSlice)); i++ {
		tempSlice := make([]*big.Int, len(tagSlice))

		for j := int64(0); j < int64(len(tagSlice)); j++ {
			var diffMonomialCoef big.Int
			diffMonomialCoef = *getDiffMonomialCoeff(tagSlice[i].XCoord, tss.EC().Params().N, j, int64(tagSlice[i].DiffTime))
			tempSlice[j] = &diffMonomialCoef
		}
		result[i] = tempSlice
	}

	if nThreshold == len(tagSlice) {
		return result, nil
	}
	tempMatrix, _ := matrix.NewMatrix(tss.EC().Params().N, result)

	resultMatrix, err := tempMatrix.DeleteColumn(uint64(nThreshold), uint64(len(tagSlice)-1))
	if err != nil {
		return emptyMatrix, err
	}
	return resultMatrix.GetMatrix(), nil

}

func GetLinearEquationCoefficient(tag Tag, degreePoly int) ([]*big.Int, error) {

	result := make([]*big.Int, degreePoly+1)

	for i := 0; i < tag.DiffTime; i++ {
		result[i] = big.NewInt(0)
	}

	for i := tag.DiffTime; i < len(result); i++ {

		result[i] = getDiffMonomialCoeff(tag.XCoord, tss.EC().Params().N, int64(i), int64(tag.DiffTime))
	}
	return result, nil
}

// Consider a monimoial x^n where n is the degree. Then ouput is n*(n_1)*...*(n-DiffTime+1)*x^{degree-DiffTimes}|_{XCoord}
// Example:x^5, DiffTime = 2 and XCoord =3 Then output is 3^(3)*5*4
func getDiffMonomialCoeff(XCoord *big.Int, fieldOrder *big.Int, degree int64, DiffTime int64) *big.Int {

	if degree < DiffTime {
		return zero
	}

	if degree == 0 {
		return one
	}

	// Get extra coefficient
	tempValue := int64(1)
	for j := int64(0); j < DiffTime; j++ {
		tempValue *= (degree - j)
	}

	var extraValue big.Int
	extraValue.SetInt64(tempValue)

	var power, result big.Int

	// x^{degree-DiffTimes}
	power.SetInt64(degree - DiffTime)
	result.Exp(XCoord, &power, fieldOrder)
	result.Mul(&result, &extraValue)
	return &result
}

// Given f(x) is a polynomial, then ouput is f^(DiffTime)(x) mod fieldOrder
// Ex: f(x)=x^5+2*x^3, DiffTime = 1, fieldOrder =3 Then f^(1)(x)= 5*x^4+6*x^2 = 2*x^4 mod 3.
func Differential(polyA Polynomial, DiffTime int) (Polynomial, error) {

	if DiffTime < 0 {
		var ZeroPoly Polynomial
		return ZeroPoly, ErrDIFFERENTIALSHOULDBEPOSITIVEINTEGER
	}

	if DiffTime >= len(polyA.coefficients) {
		var zeroPoly Polynomial
		zeroSlice := []*big.Int{big.NewInt(0)}
		zeroPoly.coefficients = zeroSlice
		return zeroPoly, nil
	}

	reduceDegree := len(polyA.coefficients) - DiffTime
	var resultPoly Polynomial
	var diffCoeffSlice = make([]*big.Int, reduceDegree)

	for i := DiffTime; i < len(polyA.coefficients); i++ {
		tempValue := int64(1)
		for j := 0; j < DiffTime; j++ {
			tempValue *= int64(i - j)
		}

		var exTra big.Int
		exTra.SetInt64(tempValue)

		var tempCoeff big.Int
		tempCoeff.Mul(polyA.coefficients[i], &exTra)
		diffCoeffSlice[i-DiffTime] = &tempCoeff
	}
	resultPoly.coefficients = diffCoeffSlice
	return resultPoly, nil
}

// Given a polynomial f(x), then the ouput is f(XCoord) mod fieldOrder
// Ex:f(3)=x^5+2*x^3, XCoord =1, fieldOrder =3 Then f(1)=3=0 mod 3
func EvaluatePolynomial(polyA Polynomial, XCoord *big.Int) *big.Int {

	if XCoord.Sign() == 0 {
		return polyA.coefficients[0]
	}

	// Compute the polynomial value using Horner's method.
	result := new(big.Int).Set(polyA.coefficients[len(polyA.coefficients)-1])
	for i := len(polyA.coefficients) - 2; i > -1; i-- {
		result.Mul(result, XCoord)
		result.Add(result, polyA.coefficients[i])
	}
	return result
}
