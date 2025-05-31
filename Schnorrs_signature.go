package main

import (
	"encoding/hex"
	"math/big"
	"strings"
)

// q = 98915E7EC8265EDFCDA31E88F24809DDB064BDC7285DD50D7289FOAC6F49DD2D
// p = EE8172AE8996608FB69359B89EB82A6985451000977A4D63BC97322CE5DC3386EA0A12B343E9190F32177539845839786BB0C345D165976EF2195EC9B1C379E3
// g = 9E96031500C8774A869582D4AFDE2127AFAD2538B4B6270A6F7C8837B50D50F206755984A49E509304D648BE2AB5AAB18EBE2CD46AC3D8495B142AA6CE23E21C

func HexToBytes(hexStr string) []byte {
	// Удаляем все пробелы и другие нежелательные символы
	clean := strings.Map(func(r rune) rune {
		switch {
		case r >= '0' && r <= '9':
			return r
		case r >= 'a' && r <= 'f':
			return r
		case r >= 'A' && r <= 'F':
			return r
		default:
			return -1 // удаляем все остальные символы
		}
	}, hexStr)

	// Декодируем в байты (паника невозможна, так как мы очистили строку)
	b, _ := hex.DecodeString(clean)
	return b
}

func BytesToBigInt(data []byte) *big.Int {
	return new(big.Int).SetBytes(data)
}

func BigIntToBytes(n *big.Int) []byte {
	return n.Bytes() // Автоматически обрабатывает nil и отрицательные числа
}

func BytesToHex(data []byte) string {
	return hex.EncodeToString(data)
}

var (
	q = BytesToBigInt(HexToBytes("98915E7EC8265EDFCDA31E88F24809DDB064BDC7285DD50D7289F0AC6F49DD2D"))
	p = BytesToBigInt(HexToBytes("EE8172AE8996608FB69359B89EB82A6985451000977A4D63BC97322CE5DC3386EA0A12B343E9190F32177539845839786BB0C345D165976EF2195EC9B1C379E3"))
	g = BytesToBigInt(HexToBytes("9E96031500C8774A869582D4AFDE2127AFAD2538B4B6270A6F7C8837B50D50F206755984A49E509304D648BE2AB5AAB18EBE2CD46AC3D8495B142AA6CE23E21C"))
)

func Schnorrs_signature(message []byte) ([]byte, []byte) {
	global := &Random_Numbers
	*global = Random_Number_Generator(*global, len(*global)+1)
	x := BytesToBigInt((*global)[len(*global)-1])
	*global = Random_Number_Generator(*global, len(*global)+1)
	r := BytesToBigInt((*global)[len(*global)-1])
	// x := BytesToBigInt(Random_Number_Generator(Random_Numbers, len(Random_Numbers))[len(Random_Numbers)])
	// fmt.Println("x: ", BytesToHex(BigIntToBytes(x)), "Length: ", len(Random_Numbers))
	// // r := BytesToBigInt(Random_Number_Generator(Random_Numbers, len(Random_Numbers))[len(Random_Numbers)])
	// fmt.Println("r: ", BytesToHex(BigIntToBytes(r)), "Length: ", len(Random_Numbers))
	P := new(big.Int).Exp(g, x, p)
	R := new(big.Int).Exp(g, r, p)
	R_byte := BigIntToBytes(R)
	P_byte := BigIntToBytes(P)
	var Base []byte
	Base = append(Base, R_byte...)
	Base = append(Base, P_byte...)
	Base = append(Base, message...)
	// fmt.Println("Base: ", BytesToHex(Hash(Base, 256)))
	E_byte := Hash(Base, 256)
	E := BytesToBigInt(E_byte)
	S := new(big.Int).Mod(new(big.Int).Add(r, new(big.Int).Mul(E, x)), q)
	S_byte := BigIntToBytes(S)
	return R_byte, S_byte
}

// func Keys_generation(message string) (*big.Int, *big.Int) {
// 	x :=
// }
