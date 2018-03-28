package main

import (
	"os"
	"fmt"
	"math"
	"math/big"
)

var (
	xLen,yLen   int64
	ctr    = 0.0
	maxctr = 0.0

	leaves    [][]*big.Int
	numslice  [][]*big.Int
	setupsl   [][]*big.Int
	solslice  [][]*big.Int
	temp      []*big.Int
	tempfloat []*big.Float
	temp1     [][]*big.Int
	bcslice   [][]*big.Int
	compslice [][]*big.Int

	//map for looking up complement to a,b,c,d splits of a number
	smap = map[float64]float64{

		0: 2,
		1: 3,
		2: 1,
		3: 0,
	}
)

func main() {

	x, y := new(big.Int), new(big.Int)

	x.SetString("3141592653589793238462643383279502884197169399375105820974944592", 10)
	y.SetString("2718281828459045235360287471352662497757247093699959574966967627", 10)

	temp = append(temp, x, y)
	numslice = append(numslice, temp)
	temp = nil

	xLen = numDigits(x)
	yLen = numDigits(y)

	if math.Ceil(math.Log2(float64(xLen))) != math.Log2(float64(xLen)) || math.Ceil(math.Log2(float64(yLen))) != math.Log2(float64(yLen)) {
		os.Exit(1)
	}

	maxctr = math.Log2(float64(xLen))

	genleaves(numslice)

	//split further if numDigits > 1

	for numDigits(setupsl[0][0]) > 1 {

		genleaves(setupsl)
		//when ctr = maxctr, the setup slice will only have numbers with 1 digit in the last leaf
		//truncating the slice to only the bottom most leaves in the tree with length 4 power #iterations
		lNum := float64(len(setupsl)) - math.Pow(4, ctr)
		setupsl = setupsl[int(lNum):]

	}

	//kick off the base case with the bottom leaves

	bcslice = baseCase(setupsl)

	//the start computing off the base case to the top node

	if bcslice != nil {

		solslice = compute(bcslice)
	} 

	//print out the final solution

	if xLen >= 64 {

		solslice[0][0].Div(solslice[0][0],big.NewInt(10))

	}


	fmt.Printf("Input Num1:%v\n", x)
	fmt.Printf("Input Num2:%v\n", y)
	fmt.Printf("Num1 & Num2 Multiplied Value by Karatsuba Recursive Algorithm:%s\n", solslice[0][0].Text(10))

	// clean up:re-initialized for next use
	setupsl = nil
	solslice = nil
	bcslice = nil

}

func genleaves(nslice [][]*big.Int) [][]*big.Int {
	var nLen, hnLen int64

	ctr++

	for s := range nslice {

		var aNum, bNum, cNum, dNum = big.NewInt(0), big.NewInt(0),
			big.NewInt(0), big.NewInt(0)
		var bigN1str, bigN2str string

		if nslice != nil {
			nLen = numDigits(nslice[0][0])
			hnLen = nLen / 2

		}

		bigN1str = nslice[s][0].Text(10)
		bigN2str = nslice[s][1].Text(10)

		aNum.SetString(bigN1str[:hnLen], 10)
		bNum.SetString(bigN1str[hnLen:], 10)
		cNum.SetString(bigN2str[:hnLen], 10)
		dNum.SetString(bigN2str[hnLen:], 10)

		//store all num splits in a temp slice until slice is full for all iterations
		//once slice is full append

		temp = append(temp, aNum, bNum, cNum, dNum)

		//set up parent slice as combinations of ac,bd,cb,da
		//uses map structure to complement pairs

		for k := range temp {
			setupsl = append(setupsl, []*big.Int{temp[k], temp[int(smap[math.Mod(float64(k), 4.0)]+math.Trunc(float64(k)/4.0)*4.0)]})
		}

		temp = nil

	}

	//return slice for base case processing
	return setupsl

}

func baseCase(leaves [][]*big.Int) [][]*big.Int {

	var end, llen int

	llen = len(leaves)
	end = int(math.Trunc(float64(llen) / 4.0))

	for j := 0; j < end; j++ {
		var acNum, bdNum, cbNum, daNum = big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(0)
		acNum = bigNmul(leaves[4*j][0], leaves[4*j][1])
		bdNum = bigNmul(leaves[4*j+1][0], leaves[4*j+1][1])
		cbNum = bigNmul(leaves[4*j+2][0], leaves[4*j+2][1])
		daNum = bigNmul(leaves[4*j+3][0], leaves[4*j+3][1])

		temp = append(temp, acNum, bdNum, cbNum, daNum)

		if len(temp) == 4 {

			bcslice = append(bcslice, temp)
			temp = nil
		}

	}

	return bcslice
}

func compute(leaves [][]*big.Int) [][]*big.Int {

	var tenPowN, tenPowNby2 = big.NewInt(0), big.NewInt(0)
	var end, llen int
	var n = new(big.Int)
	num, den := new(big.Int), new(big.Int)

	ctr--
	num.SetInt64(xLen)
	den.SetInt64(int64(math.Pow(2, ctr)))
	n.Div(num, den)
	tenPowN.Exp(big.NewInt(10), n, big.NewInt(0))
	tenPowNby2.Exp(big.NewInt(10), n.Div(n, big.NewInt(2)), big.NewInt(0))

	llen = len(leaves)
	end = int(math.Trunc(float64(llen) / 4.0))
	if end == 0 {
		end = 1
	}

	for j := 0; j < llen; j++ {
		var acNum, bdNum, cbNum, daNum, daPluscbNum, tempans = big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(0), big.NewInt(0)
		acNum = leaves[j][0]
		bdNum = leaves[j][1]
		cbNum = leaves[j][2]
		daNum = leaves[j][3]
		daPluscbNum.Add(daNum, cbNum)
		acNum = bigNmul(tenPowN, acNum)
		daPluscbNum = bigNmul(tenPowNby2, daPluscbNum)
		tempans.Add(tempans.Add(acNum, daPluscbNum), bdNum)
		temp = append(temp, tempans)

		if len(temp) == 4 || ctr == 0 {
			compslice = append(compslice, temp)
			temp = nil
		}

	}

	compslice = compslice[len(compslice)-end:]

	if ctr > 0 {
		compute(compslice)
	}

	return compslice
	

}

func bigNmul(n1 *big.Int, n2 *big.Int) *big.Int {

	var n1Pn2 = big.NewInt(0)

	n1Pn2.Mul(n1, n2)

	return n1Pn2

}

func numDigits(bigNum *big.Int) int64 {

	return int64(len(bigNum.Text(10)))

}
