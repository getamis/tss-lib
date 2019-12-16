// Copyright © 2019 AMIS Technologies
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package matrix

import (
	"math/big"

	. "github.com/onsi/ginkgo"
	. "github.com/onsi/ginkgo/extensions/table"
	. "github.com/onsi/gomega"
)

var fieldOrder = big.NewInt(101)
var bigFieldOrder, _ = new(big.Int).SetString("115792089237316195423570985008687907852837564279074904382605163141518161494337", 10)

var _ = Describe("Matrix test", func() {
	var ma = [][]*big.Int{
		{big.NewInt(1), big.NewInt(3), big.NewInt(2)},
		{big.NewInt(2), big.NewInt(7), big.NewInt(15)},
		{big.NewInt(3), big.NewInt(6), big.NewInt(81)},
	}
	var m *matrix
	BeforeEach(func() {
		var err error
		m, err = NewMatrix(fieldOrder, ma)
		Expect(err).Should(BeNil())
	})
	Context("NewMatrix()", func() {
		It("should be ok", func() {
			for i, mm := range ma {
				for j, mmm := range mm {
					Expect(m.matrix[i][j]).Should(Equal(mmm))
				}
			}
		})

		It("nil field order", func() {
			m, err := NewMatrix(nil, ma)
			Expect(err).Should(Equal(ErrZeroOrNegativeFieldOrder))
			Expect(m).Should(BeNil())
		})

		It("nil matrix", func() {
			m, err := NewMatrix(fieldOrder, nil)
			Expect(err).Should(Equal(ErrZeroRows))
			Expect(m).Should(BeNil())
		})

		It("zero column", func() {
			m, err := NewMatrix(fieldOrder, [][]*big.Int{
				{},
				{},
			})
			Expect(err).Should(Equal(ErrInconsistentColumns))
			Expect(m).Should(BeNil())
		})

		It("inconsistent column", func() {
			m, err := NewMatrix(fieldOrder, [][]*big.Int{
				{big.NewInt(1), big.NewInt(2)},
				{big.NewInt(1)},
			})
			Expect(err).Should(Equal(ErrInconsistentColumns))
			Expect(m).Should(BeNil())
		})
	})

	It("Copy", func() {
		Expect(m.Copy()).Should(Equal(m))
	})

	It("GetMatrix", func() {
		Expect(m.GetMatrix()).Should(Equal(ma))
	})

	Context("GetColumn", func() {
		It("should be ok", func() {
			for i := uint64(0); i < m.rank; i++ {
				got, err := m.GetColumn(i)
				Expect(err).Should(BeNil())

				exp := make([]*big.Int, m.rank)
				for j := uint64(0); j < m.rank; j++ {
					exp[j] = ma[j][i]
				}
				Expect(got).Should(Equal(exp))
			}
		})

		It("out of rank", func() {
			cs, err := m.GetColumn(m.rank)
			Expect(err).Should(Equal(ErrOutOfRank))
			Expect(cs).Should(BeNil())
		})
	})

	Context("GetRow", func() {
		It("should be ok", func() {
			for i := uint64(0); i < m.rank; i++ {
				got, err := m.GetRow(i)
				Expect(err).Should(BeNil())
				Expect(got).Should(Equal(ma[i]))
			}
		})

		It("out of rank", func() {
			cs, err := m.GetRow(m.rank)
			Expect(err).Should(Equal(ErrOutOfRank))
			Expect(cs).Should(BeNil())
		})
	})

	It("multiply", func() {
		matrixB, err := NewMatrix(fieldOrder, [][]*big.Int{
			{big.NewInt(1), big.NewInt(2), big.NewInt(5)},
			{big.NewInt(5), big.NewInt(8), big.NewInt(7)},
			{big.NewInt(1), big.NewInt(2), big.NewInt(5)},
		})
		Expect(err).Should(BeNil())
		got, err := m.multiply(matrixB)
		Expect(err).Should(BeNil())
		Expect(got.matrix).Should(Equal([][]*big.Int{
			{big.NewInt(18), big.NewInt(30), big.NewInt(36)},
			{big.NewInt(52), big.NewInt(90), big.NewInt(134)},
			{big.NewInt(114), big.NewInt(216), big.NewInt(462)},
		}))
	})

	Context("swapRow", func() {
		It("should be ok", func() {
			got, err := m.swapRow(0, 1)
			Expect(err).Should(BeNil())
			Expect(got.matrix).Should(Equal([][]*big.Int{
				{big.NewInt(2), big.NewInt(7), big.NewInt(15)},
				{big.NewInt(1), big.NewInt(3), big.NewInt(2)},
				{big.NewInt(3), big.NewInt(6), big.NewInt(81)},
			}))
		})

		It("out of rank", func() {
			got, err := m.swapRow(m.rank, 1)
			Expect(err).Should(Equal(ErrOutOfRank))
			Expect(got).Should(BeNil())
		})

		It("same rank", func() {
			got, err := m.swapRow(1, 1)
			Expect(err).Should(Equal(ErrSameRank))
			Expect(got).Should(BeNil())
		})
	})

	Context("swapColumn", func() {
		It("should be ok", func() {
			got, err := m.swapColumn(0, 1)
			Expect(err).Should(BeNil())
			Expect(got.matrix).Should(Equal([][]*big.Int{
				{big.NewInt(3), big.NewInt(1), big.NewInt(2)},
				{big.NewInt(7), big.NewInt(2), big.NewInt(15)},
				{big.NewInt(6), big.NewInt(3), big.NewInt(81)},
			}))
		})

		It("out of rank", func() {
			got, err := m.swapColumn(m.rank, 1)
			Expect(err).Should(Equal(ErrOutOfRank))
			Expect(got).Should(BeNil())
		})

		It("same rank", func() {
			got, err := m.swapColumn(1, 1)
			Expect(err).Should(Equal(ErrSameRank))
			Expect(got).Should(BeNil())
		})
	})

	It("Transpose", func() {
		m = m.Transpose()
		Expect(m.matrix).Should(Equal([][]*big.Int{
			{big.NewInt(1), big.NewInt(2), big.NewInt(3)},
			{big.NewInt(3), big.NewInt(7), big.NewInt(6)},
			{big.NewInt(2), big.NewInt(15), big.NewInt(81)},
		}))
	})

	DescribeTable("GetNonZeroCoefficientIndex", func(index uint64, expGot uint64, expFound bool) {
		m, err := NewMatrix(fieldOrder, [][]*big.Int{
			{big.NewInt(0), big.NewInt(18), big.NewInt(19)},
			{big.NewInt(45), big.NewInt(0), big.NewInt(81)},
			{big.NewInt(3), big.NewInt(0), big.NewInt(15)},
		})
		Expect(err).Should(BeNil())
		Expect(m).ShouldNot(BeNil())

		got, found := m.getNonZeroCoefficientIndex(index)
		Expect(got).Should(Equal(expGot))
		Expect(found).Should(Equal(expFound))
	},
		Entry("normal case", uint64(0), uint64(1), true),
		Entry("out of rank", uint64(3), uint64(0), false),
		Entry("not found", uint64(1), uint64(0), false),
	)

	It("GetGaussElimination", func() {
		gotUpper, gotLower, _, err := m.getGaussElimination()
		Expect(err).Should(BeNil())

		// Check upper matrix
		expectedUpper, err := NewMatrix(fieldOrder, [][]*big.Int{
			{big.NewInt(1), big.NewInt(3), big.NewInt(2)},
			{big.NewInt(0), big.NewInt(1), big.NewInt(11)},
			{big.NewInt(0), big.NewInt(0), big.NewInt(108)},
		})
		Expect(err).Should(BeNil())
		Expect(gotUpper.modulus().equal(expectedUpper)).Should(BeTrue())

		// Check lower matrix
		expectLower, err := NewMatrix(fieldOrder, [][]*big.Int{
			{big.NewInt(1), big.NewInt(0), big.NewInt(0)},
			{big.NewInt(-2), big.NewInt(1), big.NewInt(0)},
			{big.NewInt(-9), big.NewInt(3), big.NewInt(1)},
		})
		Expect(err).Should(BeNil())
		Expect(gotLower.modulus().equal(expectLower)).Should(BeTrue())
	})

	Context("Inverse()", func() {
		It("should be ok", func() {
			got, err := m.Inverse()
			Expect(err).Should(BeNil())
			Expect(got.GetMatrix()).Should(Equal([][]*big.Int{
				{big.NewInt(97), big.NewInt(68), big.NewInt(91)},
				{big.NewInt(41), big.NewInt(54), big.NewInt(85)},
				{big.NewInt(42), big.NewInt(87), big.NewInt(29)},
			}))
		})

		It("should be ok with big prime", func() {
			m, err := NewMatrix(bigFieldOrder, ma)
			Expect(err).Should(BeNil())
			Expect(m).ShouldNot(BeNil())

			got, err := m.Inverse()
			Expect(err).Should(BeNil())

			expected := make([][]*big.Int, 3)
			for i := 0; i < 3; i++ {
				expected[i] = make([]*big.Int, 3)
			}
			expected[0][0], _ = new(big.Int).SetString("67545385388434447330416407921734612914155245829460360889853011832552260871701", 10)
			expected[0][1], _ = new(big.Int).SetString("93276960774504712980098849034776370214785815669254784085987492530667407870436", 10)
			expected[0][2], _ = new(big.Int).SetString("31092320258168237660032949678258790071595271889751594695329164176889135956813", 10)
			expected[1][0], _ = new(big.Int).SetString("9649340769776349618630915417390658987736463689922908698550430261793180124527", 10)
			expected[1][1], _ = new(big.Int).SetString("112575642314057412217360679869557688190258743049100601483088353054253768119495", 10)
			expected[1][2], _ = new(big.Int).SetString("37525214104685804072453559956519229396752914349700200494362784351417922706498", 10)
			expected[2][0], _ = new(big.Int).SetString("9649340769776349618630915417390658987736463689922908698550430261793180124528", 10)
			expected[2][1], _ = new(big.Int).SetString("73978279234952013742837018199995052239312888289408966688886632007081047621382", 10)
			expected[2][2], _ = new(big.Int).SetString("101854152569861468196659662739123622648329338949186258484698986096705790203352", 10)
			Expect(got.matrix).Should(Equal(expected))
		})

		It("should be ok with big prime 2", func() {
			m, err := NewMatrix(bigFieldOrder, [][]*big.Int{
				{big.NewInt(4), big.NewInt(10), big.NewInt(30)},
				{big.NewInt(10), big.NewInt(30), big.NewInt(100)},
				{big.NewInt(30), big.NewInt(100), big.NewInt(354)},
			})
			Expect(err).Should(BeNil())
			Expect(m).ShouldNot(BeNil())

			got, err := m.Inverse()
			Expect(err).Should(BeNil())

			expected := make([][]*big.Int, 3)
			for i := 0; i < 3; i++ {
				expected[i] = make([]*big.Int, 3)
			}
			expected[0][0], _ = new(big.Int).SetString("28948022309329048855892746252171976963209391069768726095651290785379540373592", 10)
			expected[0][1], _ = new(big.Int).SetString("86844066927987146567678238756515930889628173209306178286953872356138621120746", 10)
			expected[0][2], _ = new(big.Int).SetString("86844066927987146567678238756515930889628173209306178286953872356138621120754", 10)
			expected[1][0], _ = new(big.Int).SetString("86844066927987146567678238756515930889628173209306178286953872356138621120746", 10)
			expected[1][1], _ = new(big.Int).SetString("17368813385597429313535647751303186177925634641861235657390774471227724224157", 10)
			expected[1][2], _ = new(big.Int).SetString("28948022309329048855892746252171976963209391069768726095651290785379540373583", 10)
			expected[2][0], _ = new(big.Int).SetString("86844066927987146567678238756515930889628173209306178286953872356138621120754", 10)
			expected[2][1], _ = new(big.Int).SetString("28948022309329048855892746252171976963209391069768726095651290785379540373583", 10)
			expected[2][2], _ = new(big.Int).SetString("86844066927987146567678238756515930889628173209306178286953872356138621120753", 10)
			Expect(got.matrix).Should(Equal(expected))
		})

		It("inverse not exist", func() {
			m, err := NewMatrix(fieldOrder, [][]*big.Int{
				{big.NewInt(0), big.NewInt(3), big.NewInt(2)},
				{big.NewInt(0), big.NewInt(7), big.NewInt(15)},
				{big.NewInt(0), big.NewInt(6), big.NewInt(81)},
			})
			Expect(err).Should(BeNil())
			Expect(m).ShouldNot(BeNil())
			got, err := m.Inverse()
			Expect(err).Should(Equal(ErrNotInvertableMatrix))
			Expect(got).Should(BeNil())
		})
	})

	It("multiInverseDiagonal", func() {
		m, err := NewMatrix(bigFieldOrder, [][]*big.Int{
			{big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(5)},
			{big.NewInt(2), big.NewInt(0), big.NewInt(0), big.NewInt(0)},
			{big.NewInt(0), big.NewInt(3), big.NewInt(0), big.NewInt(0)},
			{big.NewInt(0), big.NewInt(0), big.NewInt(-1), big.NewInt(0)},
		})
		Expect(err).Should(BeNil())
		Expect(m).ShouldNot(BeNil())

		diagonalMatrice, _ := NewMatrix(bigFieldOrder, [][]*big.Int{
			{big.NewInt(2), big.NewInt(0), big.NewInt(0), big.NewInt(0)},
			{big.NewInt(0), big.NewInt(3), big.NewInt(0), big.NewInt(0)},
			{big.NewInt(0), big.NewInt(0), big.NewInt(-1), big.NewInt(0)},
			{big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(14)},
		})

		got, err := m.multiInverseDiagonal(diagonalMatrice)
		Expect(err).Should(BeNil())

		inverse11 := new(big.Int).ModInverse(big.NewInt(2), bigFieldOrder)
		inverse22 := new(big.Int).ModInverse(big.NewInt(3), bigFieldOrder)
		inverse33 := new(big.Int).ModInverse(big.NewInt(-1), bigFieldOrder)
		inverse44 := new(big.Int).ModInverse(big.NewInt(14), bigFieldOrder)

		original, err := NewMatrix(bigFieldOrder, [][]*big.Int{
			{big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(5)},
			{big.NewInt(2), big.NewInt(0), big.NewInt(0), big.NewInt(0)},
			{big.NewInt(0), big.NewInt(3), big.NewInt(0), big.NewInt(0)},
			{big.NewInt(0), big.NewInt(0), big.NewInt(-1), big.NewInt(0)},
		})
		Expect(err).Should(BeNil())

		inverseDiagonal, err := NewMatrix(bigFieldOrder, [][]*big.Int{
			{inverse11, big.NewInt(0), big.NewInt(0), big.NewInt(0)},
			{big.NewInt(0), inverse22, big.NewInt(0), big.NewInt(0)},
			{big.NewInt(0), big.NewInt(0), inverse33, big.NewInt(0)},
			{big.NewInt(0), big.NewInt(0), big.NewInt(0), inverse44},
		})
		Expect(err).Should(BeNil())

		expected, err := inverseDiagonal.multiply(original)
		Expect(err).Should(BeNil())
		Expect(got.equal(expected)).Should(BeTrue())
	})

	It("modInverse", func() {
		m, err := NewMatrix(bigFieldOrder, [][]*big.Int{
			{big.NewInt(0), big.NewInt(1), big.NewInt(4)},
			{big.NewInt(0), big.NewInt(3), big.NewInt(5)},
			{big.NewInt(0), big.NewInt(0), big.NewInt(2)},
		})
		Expect(err).Should(BeNil())
		Expect(m).ShouldNot(BeNil())
		got := m.modInverse(1, 2)
		expected := new(big.Int).ModInverse(big.NewInt(5), bigFieldOrder)
		Expect(got.Cmp(expected)).Should(BeZero())
	})
})
