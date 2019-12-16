package matrix

import (
	"errors"
	"math/big"
)

var (
	// big0 is the constant 0
	big0 = big.NewInt(0)

	// ErrZeroOrNegativeFieldOrder is returned if the field order is zero or negative
	ErrZeroOrNegativeFieldOrder = errors.New("zero or negative field order")
	// ErrNilMatrix is returned if it's a nil matrix
	ErrNilMatrix = errors.New("nil matrix")
	// ErrZeroRows is returned if the number of row of the matrix is zero
	ErrZeroRows = errors.New("zero rows")
	// ErrZeroColumns is returned if the number of column of the matrix is zero
	ErrZeroColumns = errors.New("zero columns")
	// ErrInconsistentColumns is returned if the column rank is inconsistent in this matrix
	ErrInconsistentColumns = errors.New("inconsistent columns")
	// ErrZeroOrNegativeRank is returned if the rank is zero or negative
	ErrZeroOrNegativeRank = errors.New("zero or negative rank")
	// ErrOutOfRank is returned if the index is out of the column or row rank
	ErrOutOfRank = errors.New("out of rank")
	// ErrSameRank is returned if the ranks are the same
	ErrSameRank = errors.New("same rank")
	// ErrInconsistentRank is returned if the two matrixes are the inconsistent rank
	ErrInconsistentRank = errors.New("inconsistent rank")
	// ErrNotSquareMatrix is returned if it's not a square matrix
	ErrNotSquareMatrix = errors.New("not a square matrix")
	// ErrNotInvertableMatrix is returned if it's not an invertable matrix
	ErrNotInvertableMatrix = errors.New("not invertable matrix")
	// ErrRangeIndex is returned if any index is false
	ErrRangeIndex = errors.New("The range of index is false")
)

// Matrix is the struct for matrix operation
type Matrix struct {
	fieldOrder *big.Int

	// number of row of this matrix
	numberRow uint64

	// number of column of this matrix
	numberColumn uint64
	matrix       [][]*big.Int
}

// NewMatrix checks the input matrix slices. It returns error if the
// rank of rows or columns is zero or the rank of column is inconsistent.
func NewMatrix(fieldOrder *big.Int, matrix [][]*big.Int) (*Matrix, error) {
	if fieldOrder == nil || fieldOrder.Cmp(big.NewInt(0)) <= 0 {
		return nil, ErrZeroOrNegativeFieldOrder
	}
	numberRow := uint64(len(matrix))
	if numberRow == 0 {
		return nil, ErrZeroRows
	}
	numberColumn := uint64(len(matrix[0]))
	if numberColumn == 0 {
		return nil, ErrZeroColumns
	}
	for i := uint64(0); i < numberRow; i++ {
		if uint64(len(matrix[i])) != numberColumn {
			return nil, ErrInconsistentColumns
		}
		for j := uint64(0); j < numberColumn; j++ {
			if matrix[i][j] == nil {
				return nil, ErrNilMatrix
			}
		}
	}
	return &Matrix{
		fieldOrder:   fieldOrder,
		numberRow:    numberRow,
		numberColumn: numberColumn,
		matrix:       matrix,
	}, nil
}

// Copy returns a copied matrix
func (m *Matrix) Copy() *Matrix {
	return &Matrix{
		fieldOrder:   new(big.Int).Set(m.fieldOrder),
		numberRow:    m.numberRow,
		numberColumn: m.numberColumn,
		matrix:       m.GetMatrix(),
	}
}

func (m *Matrix) GetMatrix() [][]*big.Int {

	newMatrix := make([][]*big.Int, m.numberRow)
	for i := uint64(0); i < m.numberRow; i++ {
		newMatrix[i] = make([]*big.Int, m.numberColumn)
		for j := uint64(0); j < m.numberColumn; j++ {
			newMatrix[i][j] = m.Get(i, j)
		}
	}
	return m.matrix
}

// GetColumn gets column at the index
// Assume matrixA = [ 1,2,3 ]
//                  [ 2,4,5 ]
//                  [ 5,10,3]
// Then the output of GetColumn(matrixA, nIndex ) is the indicated column.
// Ex: GetColumn(matrixA, 2 )= [3,5,3], GetColumn(matrixA, 1 )=[2,4,10]
func (m *Matrix) GetColumn(nIndex uint64) ([]*big.Int, error) {
	if nIndex >= m.numberColumn {
		return nil, ErrOutOfRank
	}

	tempSlice := make([]*big.Int, m.numberRow)
	for i := uint64(0); i < m.numberRow; i++ {
		tempSlice[i] = m.Get(i, nIndex)
	}
	return tempSlice, nil
}

// GetRow gets row at the index
// Assume matrixA = [ 1,2,3 ]
//                  [ 2,4,5 ]
//                  [ 5,10,3]
// Then the ouput of GetColumn(matrixA, nIndex ) is the indicated row.
// Ex: GetRow(matrixA, 2 )= [5,10,3], GetRow(matrixA, 1 )=[2,4,5]
func (m *Matrix) GetRow(nIndex uint64) ([]*big.Int, error) {
	if nIndex >= m.numberRow {
		return nil, ErrOutOfRank
	}
	tempSlice := make([]*big.Int, m.numberColumn)
	for i := uint64(0); i < m.numberColumn; i++ {
		tempSlice[i] = m.Get(nIndex, i)
	}
	return tempSlice, nil
}

// Get gets the element at (i, j)
func (m *Matrix) Get(i, j uint64) *big.Int {
	v := m.get(i, j)
	if v == nil {
		return nil
	}
	return new(big.Int).Mod(v, m.fieldOrder)
}

// get gets the element at (i, j) without mod its value
func (m *Matrix) get(i, j uint64) *big.Int {
	if i >= m.numberRow {
		return nil
	}
	if j >= m.numberColumn {
		return nil
	}
	return m.matrix[i][j]
}

// get gets the element at (i, j) without checking the index range and ensure the returned value is in our field
func (m *Matrix) modInverse(i, j uint64) *big.Int {
	v := m.get(i, j)
	if v == nil {
		return nil
	}
	return new(big.Int).ModInverse(v, m.fieldOrder)
}

// Transpose transposes the matrix
// This function give the transpose of input.
// Ex: A =[ 1, 2  ] (i.e. 1X2 matrix)
// output is [ 1 ] (i.e. 2X1 matrix)
//          [ 2 ]
func (m *Matrix) Transpose() *Matrix {
	transposeMatrix := make([][]*big.Int, m.numberColumn)
	for i := uint64(0); i < m.numberColumn; i++ {
		tempSlice := make([]*big.Int, m.numberRow)

		for j := uint64(0); j < m.numberRow; j++ {
			tempSlice[j] = m.matrix[j][i]
		}
		transposeMatrix[i] = tempSlice
	}
	m.matrix = transposeMatrix

	// Exchange rank
	tempRank := m.numberColumn
	m.numberColumn = m.numberRow
	m.numberRow = tempRank
	return m
}

// Add adds the matrix
// The standard addition of Matrices
func (m *Matrix) Add(matrix *Matrix) (*Matrix, error) {
	if matrix == nil {
		return nil, ErrNilMatrix
	}
	if m.numberColumn != matrix.numberColumn || m.numberRow != matrix.numberRow {
		return nil, ErrInconsistentRank
	}

	for i := uint64(0); i < m.numberRow; i++ {
		m.matrix[i] = addSlices(m.matrix[i], matrix.matrix[i])
	}
	return m.modulus(), nil
}

func (m *Matrix) multiply(matrix *Matrix) (*Matrix, error) {
	if matrix == nil {
		return nil, ErrNilMatrix
	}
	// check two matrices can do multiplication by checking their rank
	if m.numberColumn != matrix.numberRow {
		return nil, ErrInconsistentRank
	}

	for i := uint64(0); i < m.numberRow; i++ {
		tempSlice := make([]*big.Int, matrix.numberColumn)
		for j := uint64(0); j < matrix.numberColumn; j++ {
			tempValue := big.NewInt(0)
			for k := uint64(0); k < m.numberColumn; k++ {
				tempValue.Add(tempValue, new(big.Int).Mul(m.matrix[i][k], matrix.matrix[k][j]))
			}
			tempSlice[j] = tempValue
		}
		m.matrix[i] = tempSlice
	}
	m.numberColumn = matrix.numberColumn
	return m, nil
}

// All components of a matrix modulus a fieldOrder.
// Ex: A = [10,9]    and fieldOrder = 7
//         [23,14]
// Then output is [3,2]
//               [2,0]
func (m *Matrix) modulus() *Matrix {
	for i := uint64(0); i < m.numberRow; i++ {
		for j := uint64(0); j < m.numberColumn; j++ {
			m.matrix[i][j].Mod(m.matrix[i][j], m.fieldOrder)
		}
	}
	return m
}

// Interchange two rows of a given matrix.
// Ex: A = [10,9]    and fieldOrder = 7
//         [23,14]
// SwapRow(A,0,1) = [23,14]
//                  [10,9]
func (m *Matrix) swapRow(nIndexRow1 uint64, nIndexRow2 uint64) (*Matrix, error) {
	if m.numberRow < nIndexRow1 || m.numberRow < nIndexRow2 {
		return nil, ErrOutOfRank
	}

	for i := uint64(0); i < m.numberColumn; i++ {
		m.matrix[nIndexRow1][i], m.matrix[nIndexRow2][i] = m.matrix[nIndexRow2][i], m.matrix[nIndexRow1][i]
	}
	return m, nil
}

func (m *Matrix) swapColumn(nIndexColumn1 uint64, nIndexColumn2 uint64) (*Matrix, error) {
	if m.numberColumn < nIndexColumn1 || m.numberColumn < nIndexColumn2 {
		return nil, ErrOutOfRank
	}

	for i := uint64(0); i < m.numberRow; i++ {
		m.matrix[i][nIndexColumn1], m.matrix[i][nIndexColumn2] = m.matrix[i][nIndexColumn2], m.matrix[i][nIndexColumn1]
	}
	return m, nil
}

// If this matrix is square then ouput its rank; Otherwise, ouput errors.
func (m *Matrix) IsSquare() bool {
	return m.numberColumn == m.numberRow
}

// Inverse Matrix
// Get the inverse matrix of A
func (m *Matrix) Inverse() (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquareMatrix
	}
	// Get U, L^{-1}. Note that A= L*U
	upperMatrix, lowerMatrix, _, err := m.getGaussElimination()

	if err != nil {
		return nil, err
	}

	copyLowerMatrix := lowerMatrix.Copy()
	// K=U^t
	upperMatrix.Transpose()
	// Get D, L_K^{-1}. Note that K=L_K*D
	tempUpperResult, tempLowerResult, _, err := upperMatrix.getGaussElimination()
	if err != nil {
		return nil, err
	}

	tempResult, err := tempLowerResult.multiInverseDiagonal(tempUpperResult)
	if err != nil {
		return nil, err
	}
	// Get (D^{-1}L_{K}^{-1})^t = ((L_K*D)^{-1})^t = (K^{-1})^{t}, so the transpose of (K^{-1})^{t} is U^{-1}
	tempResult.Transpose()

	// U^{-1}*L^{-1} = (L*U)^{-1} = A^{-1}
	tempResult, err = tempResult.multiply(copyLowerMatrix)
	if err != nil {
		return nil, err
	}
	m = tempResult.modulus()
	return m, nil
}

// Determinant returns the determinant of the matrix
func (m *Matrix) Determinant() (*big.Int, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquareMatrix
	}

	m.modulus()
	// We only use elementary matrix (i.e. its determine is 1), so det(upperMatrix)=det(A).
	// Furthermore, upperMatrix is a uppertriangular matrix. Thus, the determinant of this matrix
	// is the multiplication of all digonal elements.
	upperMatrix, _, permutationTimes, err := m.getGaussElimination()
	if err != nil {
		return big.NewInt(0), nil
	}

	rank := m.numberRow
	result := big.NewInt(1)
	for i := uint64(0); i < rank; i++ {
		result.Mul(result, upperMatrix.matrix[i][i])
	}

	// negative result if the times of permutation is odd
	if permutationTimes%2 == 1 {
		result.Neg(result)
	}

	result.Mod(result, m.fieldOrder)
	return result, nil
}

// As give the index of row of a matrix, this function will find nonzero value such that this value has the smallest index of row.
func (m *Matrix) getNonZeroCoefficientIndex(fromRowIndex uint64) (uint64, bool) {
	for i := fromRowIndex; i < m.numberRow; i++ {
		if m.Get(i, fromRowIndex).Cmp(big.NewInt(0)) != 0 {
			return i, true
		}
	}
	return 0, false
}

// Only work "matrixA is squared-matrix"
// Then the output is U_A and L^{-1} such that LU_A = A. Here U_A is a upper triangular matrix
// with det(U_A) = det(A). (i.e. <A|I> = <U_A|L^{-1}> by Gauss elimination)
func (m *Matrix) getGaussElimination() (*Matrix, *Matrix, int, error) {
	if !m.IsSquare() {
		return nil, nil, 0, ErrNotSquareMatrix
	}
	rank := m.numberRow
	lower, err := newIdentityMatrix(rank, m.fieldOrder)
	if err != nil {
		return nil, nil, 0, err
	}

	upper := m.Copy()
	permutationTimes := 0

	for i := uint64(0); i < rank; i++ {
		changeIndex, found := upper.getNonZeroCoefficientIndex(i)
		if !found {
			return nil, nil, 0, ErrNotInvertableMatrix
		}

		// If the index is changed, swap rows
		if i != changeIndex {
			permutationTimes++
			// Swap lower and higher matrix
			upper, err = upper.swapRow(i, changeIndex)
			if err != nil {
				return nil, nil, 0, err
			}

			lower, err = lower.swapRow(i, changeIndex)
			if err != nil {
				return nil, nil, 0, err
			}
		}

		inverse := upper.modInverse(i, i)
		if inverse == nil {
			return nil, nil, 0, ErrNotInvertableMatrix
		}
		for j := i + 1; j < rank; j++ {
			tempValue := new(big.Int).Mul(upper.matrix[j][i], inverse)
			inverseDiagonalComponent := new(big.Int).Neg(tempValue)

			// Make (j, i) element to zero at upper matrix
			rowI, err := upper.GetRow(i)
			if err != nil {
				return nil, nil, 0, err
			}
			rowJ, err := upper.GetRow(j)
			if err != nil {
				return nil, nil, 0, err
			}
			tempResultASlice := multiScalar(rowI, inverseDiagonalComponent)
			upper.matrix[j] = addSlices(rowJ, tempResultASlice)

			// Do the same above operation at lower matrix
			rowLowerI, err := lower.GetRow(i)
			if err != nil {
				return nil, nil, 0, err
			}
			rowLowerJ, err := lower.GetRow(j)
			if err != nil {
				return nil, nil, 0, err
			}

			tempResultIdentitySlice := multiScalar(rowLowerI, inverseDiagonalComponent)
			lower.matrix[j] = addSlices(rowLowerJ, tempResultIdentitySlice)
		}
	}
	// TODO: We can remove it.
	upper = upper.modulus()
	lower = lower.modulus()
	return upper, lower, permutationTimes, nil
}

func (m *Matrix) exchangeRowOrColumnNonZeroElement(indexDiagonal uint64) (*Matrix, bool) {
	for i := indexDiagonal; i < m.numberRow; i++ {
		if m.matrix[i][indexDiagonal].Cmp(big.NewInt(0)) != 0 {
			m.swapRow(i, indexDiagonal)
			return m, true
		}
	}

	for i := indexDiagonal; i < m.numberColumn; i++ {
		if m.matrix[indexDiagonal][i].Cmp(big.NewInt(0)) != 0 {
			m.swapColumn(i, indexDiagonal)
			return m, true
		}
	}
	return m, false
}

// if fieldOrder <= 0 then return the standard rank
// Otherwise, return the rank in the ring Z/fieldOrderZ
func (m *Matrix) GetMatrixRank(fieldOrder *big.Int) (uint64, error) {
	upper := m.Copy()
	if upper.numberRow < upper.numberColumn {
		upper = upper.Transpose()
	}

	for i := uint64(0); i < upper.numberColumn; i++ {
		upper, found := upper.exchangeRowOrColumnNonZeroElement(i)
		if !found {
			continue
		}

		inverse := upper.modInverse(i, i)

		if inverse == nil {
			continue
		}

		for j := i + 1; j < upper.numberRow; j++ {
			tempValue := new(big.Int).Mul(upper.matrix[j][i], inverse)
			inverseDiagonalComponent := new(big.Int).Neg(tempValue)

			// Make (j, i) element to zero at upper matrix
			rowI, err := upper.GetRow(i)
			if err != nil {
				return 0, err
			}
			rowJ, err := upper.GetRow(j)
			if err != nil {
				return 0, err
			}
			tempResultASlice := multiScalar(rowI, inverseDiagonalComponent)
			upper.matrix[j] = addSlices(rowJ, tempResultASlice)
		}
	}

	result := uint64(0)
	// Compute rank
	if fieldOrder.Cmp(big.NewInt(0)) < 1 {
		for i := uint64(0); i < upper.numberColumn; i++ {
			if upper.matrix[i][i].Cmp(big.NewInt(0)) != 0 {
				result++
			}
		}
		return result, nil
	}
	for i := uint64(0); i < upper.numberColumn; i++ {
		var tempBig big.Int
		tempBig.Mod(upper.matrix[i][i], fieldOrder)
		if tempBig.Cmp(big.NewInt(0)) != 0 {
			result++
		}
	}
	return result, nil
}

// multiInverseDiagonal inverse the diagonal matrix and multiplies it
// Only use in computing inverse matrix.
func (m *Matrix) multiInverseDiagonal(diagonal *Matrix) (*Matrix, error) {
	rank := m.numberRow
	for i := uint64(0); i < rank; i++ {
		inverse := diagonal.modInverse(i, i)
		if inverse == nil {
			return nil, ErrNotInvertableMatrix
		}
		for j := uint64(0); j < rank; j++ {
			m.matrix[i][j].Mul(m.matrix[i][j], inverse)
		}
	}
	return m, nil
}

// We will delete the row from nLowerIndex to nUpperIndex
// Ex:
// a_11 a_12 a_13 a_14
// a_21 a_22 a_23 a_24
// a_31 a_32 a_33 a_34
// a_41 a_42 a_43 a_44
// Then DeleteRow(1,2) will gives
// a_11 a_12 a_13 a_14
// a_41 a_42 a_43 a_44
func (m *Matrix) DeleteRow(nLowerIndex, nUpperIndex uint64) (*Matrix, error) {
	var emptyMatrix *Matrix
	if nLowerIndex < 0 || nUpperIndex+1 > m.numberRow {
		return emptyMatrix, ErrRangeIndex
	}

	var reduceMatrix [][]*big.Int

	for i := uint64(0); i < nLowerIndex; i++ {
		reduceMatrix = append(reduceMatrix, m.matrix[i])
	}

	for i := nUpperIndex + 1; i < m.numberRow; i++ {
		reduceMatrix = append(reduceMatrix, m.matrix[i])
	}
	resultMatrix, err := NewMatrix(m.fieldOrder, reduceMatrix)

	if err != nil {
		return emptyMatrix, err
	}
	return resultMatrix, nil
}

func (m *Matrix) DeleteColumn(nLowerIndex, nUpperIndex uint64) (*Matrix, error) {
	var emptyMatrix *Matrix
	transposeM := m.Transpose()

	if nLowerIndex < 0 || nUpperIndex+1 > transposeM.numberRow {
		return emptyMatrix, ErrRangeIndex
	}

	var reduceMatrix [][]*big.Int

	for i := uint64(0); i < nLowerIndex; i++ {
		reduceMatrix = append(reduceMatrix, transposeM.matrix[i])
	}

	for i := nUpperIndex + 1; i < transposeM.numberRow; i++ {
		reduceMatrix = append(reduceMatrix, transposeM.matrix[i])
	}
	resultMatrix, err := NewMatrix(transposeM.fieldOrder, reduceMatrix)

	if err != nil {
		return emptyMatrix, err
	}
	resultMatrix.Transpose()
	return resultMatrix, nil
}

// This is a special case of Pseudoinverse
// If the columns of m are linearly independent (so that row rank >= column rank)
// If m^t*m is invertible. In this case, an explicity formula is : (m^t*m)^(-1)*m^t.
func (m *Matrix) Pseudoinverse() (*Matrix, error) {
	copy := m.Copy()
	copyTranspose := m.Copy()
	copyTranspose.Transpose()

	copyTran := m.Copy()
	copyTran.Transpose()

	symmetricForm, err := copyTranspose.multiply(copy)
	if err != nil {
		return m, err
	}

	// (m^t*m)^(-1)
	inverseSymmetric, err := symmetricForm.Inverse()
	if err != nil {
		return m, err
	}

	// (m^t*m)^(-1)*m^t
	result, err := inverseSymmetric.multiply(copyTran)
	result.modulus()

	if err != nil {
		return m, err
	}

	return result, nil
}

func (m *Matrix) Equal(m2 *Matrix) bool {
	if m2 == m {
		return true
	}

	if m.numberRow != m2.numberRow {
		return false
	}

	if m.numberColumn != m2.numberColumn {
		return false
	}
	if m.fieldOrder.Cmp(m2.fieldOrder) != 0 {
		return false
	}

	for i, mm := range m.matrix {
		for j := range mm {
			if m.Get(uint64(i), uint64(j)).Cmp(m2.Get(uint64(i), uint64(j))) != 0 {
				return false
			}
		}
	}

	return true
}
