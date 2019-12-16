// Copyright © 2019 Binance
//
// This file is part of Binance. The full Binance copyright notice, including
// terms governing use, modification, and redistribution, is contained in the
// file LICENSE at the root of the source code distribution tree.

package signing

import (
	"fmt"
	"math/big"

	"github.com/binance-chain/tss-lib/common"
	"github.com/binance-chain/tss-lib/crypto"
	"github.com/binance-chain/tss-lib/crypto/vss"
	"github.com/binance-chain/tss-lib/tss"
)

// PrepareForSigning(), GG18Spec (11) Fig. 14
func PrepareForSigning(i, pax, threshold int, xi *big.Int, ks []*big.Int, ranks []uint, bigXs []*crypto.ECPoint) (wi *big.Int, bigWs []*crypto.ECPoint) {
	modQ := common.ModInt(tss.EC().Params().N)
	if len(ks) != len(bigXs) {
		panic(fmt.Errorf("indices and bigX are not same length"))
	}
	if len(ks) != pax {
		panic(fmt.Errorf("indices is not in pax size"))
	}

	tagSlice := make([]vss.Tag, len(ks))

	for j := 0; j < len(tagSlice); j++ {
		tagSlice[j].XCoord = ks[j]
		tagSlice[j].DiffTime = int(ranks[j])
	}
	birkhoff, _ := vss.GetBirkhoffCiefficient(tagSlice, threshold + 1)

	wi = xi
	wi = new(big.Int).Mul(wi, birkhoff[i])
	wi.Mod(wi, tss.EC().Params().N)

	//fmt.Println(wi)
	// 2-4.

	wi = xi
	for j := 0; j < pax; j++ {
		if j == i {
			continue
		}
		// big.Int Div is calculated as: a/b = a * modInv(b,q)
		coef := modQ.Mul(ks[j], modQ.ModInverse(new(big.Int).Sub(ks[j], ks[i])))
		wi = modQ.Mul(wi, coef)
	}
	//fmt.Println(wi)

	// 5-10.
	bigWs = make([]*crypto.ECPoint, len(ks))
	for j := 0; j < pax; j++ {
		bigWj := bigXs[j]
		bigWj = bigWj.ScalarMult(birkhoff[j])
		/*
			for c := 0; c < pax; c++ {
				if j == c {
					continue
				}
				ksc := ks[c]
				ksj := ks[j]
				if ksj.Cmp(ksc) == 0 {
					panic(fmt.Errorf("index of two parties are equal"))
				}
				// big.Int Div is calculated as: a/b = a * modInv(b,q)
				iota := modQ.Mul(ksc, modQ.ModInverse(new(big.Int).Sub(ksc, ksj)))
				bigWj = bigWj.ScalarMult(iota)
			}
		*/
		bigWs[j] = bigWj
	}
	return
}
