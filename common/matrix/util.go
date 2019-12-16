package matrix

import (
	"math/big"
)

// Generate the identity matrix with rank nRank
func newIdentityMatrix(rank uint64, fieldOrder *big.Int) (*Matrix, error) {
	if rank < 1 {
		return nil, ErrZeroOrNegativeRank
	}
	identityMatrix := make([][]*big.Int, rank)

	for i := uint64(0); i < rank; i++ {
		tempSlice := make([]*big.Int, rank)

		for j := uint64(0); j < rank; j++ {
			tempSlice[j] = big.NewInt(0)
		}

		tempSlice[i] = big.NewInt(1)
		identityMatrix[i] = tempSlice
	}
	return &Matrix{
		fieldOrder:   fieldOrder,
		numberRow:    rank,
		numberColumn: rank,
		matrix:       identityMatrix,
	}, nil
}

// MultiScalar multiply the slice with the scalar
// We assume that the len(slice) != 0
func multiScalar(slice []*big.Int, scalar *big.Int) []*big.Int {
	result := make([]*big.Int, len(slice))
	for i := 0; i < len(slice); i++ {
		result[i] = new(big.Int).Mul(slice[i], scalar)
	}
	return result
}

// AddSlices adds two slices
// We assume that the len(sliceA) = len(sliceB) != 0
func addSlices(sliceA []*big.Int, sliceB []*big.Int) []*big.Int {
	result := make([]*big.Int, len(sliceA))
	for i := 0; i < len(sliceA); i++ {
		result[i] = new(big.Int).Add(sliceA[i], sliceB[i])
	}
	return result
}
