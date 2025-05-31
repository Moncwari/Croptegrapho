package main

import (
	// "bufio"
	"encoding/binary"
	"encoding/hex"

	// "fmt"
	"math/big"
	// "os"
	// "strings"
)

func Vec_n(z int, n int) []byte {
	vec := make([]byte, n)
	for i := 0; i < n; i++ {
		vec[i] = byte(z & 1) // Получаем младший бит
		z >>= 1              // Сдвиг вправо (деление на 2)
	}
	return vec
}

// Int_n преобразует бинарный вектор длины n обратно в целое число
func Int_n(vec []byte) int {
	z := 0
	for i := len(vec) - 1; i >= 0; i-- {
		z = (z << 1) | int(vec[i]) // Сдвиг влево (умножение на 2) и добавление бита
	}
	return z
}

func getBit(b byte, pos uint) byte {
	return (b >> pos) & 1
}

// L performs a linear transformation on a 512-bit input block V.
// The transformation involves multiplying each 64-bit segment of V
// by a predefined matrix and XORing the results. The function returns
// a new 512-bit block resulting from the transformation.
func L(V []byte) []byte {
	var Byte_matrix = [64]uint64{
		0x8e20faa72ba0b470, 0x47107ddd9b505a38, 0xad08b0e0c3282d1c, 0xd8045870ef14980e,
		0x6c022c38f90a4c07, 0x3601161cf205268d, 0x1b8e0b0e798c13c8, 0x83478b07b2468764,
		0xa011d380818e8f40, 0x5086e740ce47c920, 0x2843fd2067adea10, 0x14aff010bdd87508,
		0x0ad97808d06cb404, 0x05e23c0468365a02, 0x8c711e02341b2d01, 0x46b60f011a83988e,
		0x90dab52a387ae76f, 0x486dd4151c3dfdb9, 0x24b86a840e90f0d2, 0x125c354207487869,
		0x092e94218d243cba, 0x8a174a9ec8121e5d, 0x4585254f64090fa0, 0xaccc9ca9328a8950,
		0x9d4df05d5f661451, 0xc0a878a0a1330aa6, 0x60543c50de970553, 0x302a1e286fc58ca7,
		0x18150f14b9ec46dd, 0x0c84890ad27623e0, 0x0642ca05693b9f70, 0x0321658cba93c138,
		0x86275df09ce8aaa8, 0x439da0784e745554, 0xafc0503c273aa42a, 0xd960281e9d1d5215,
		0xe230140fc0802984, 0x71180a8960409a42, 0xb60c05ca30204d21, 0x5b068c651810a89e,
		0x456c34887a3805b9, 0xac361a443d1c8cd2, 0x561b0d22900e4669, 0x2b838811480723ba,
		0x9bcf4486248d9f5d, 0xc3e9224312c8c1a0, 0xeffa11af0964ee50, 0xf97d86d98a327728,
		0xe4fa2054a80b329c, 0x727d102a548b194e, 0x39b008152acb8227, 0x9258048415eb419d,
		0x492c024284fbaec0, 0xaa16012142f35760, 0x550b8e9e21f7a530, 0xa48b474f9ef5dc18,
		0x70a6a56e2440598e, 0x3853dc371220a247, 0x1ca76e95091051ad, 0x0edd37c48a08a6d8,
		0x07e095624504536c, 0x8d70c431ac02a736, 0xc83862965601dd1b, 0x641c314b2b8ee083,
	}

	L := make([]byte, 64)

	for l := 0; l < 8; l++ {
		b := binary.BigEndian.Uint64(V[l*8 : (l+1)*8])
		var c uint64 = 0

		for i := 0; i < 64; i++ {
			if (b>>(63-i))&1 == 1 {
				c ^= Byte_matrix[i]
			}
		}

		// Записываем результат обратно в `L`
		binary.BigEndian.PutUint64(L[l*8:(l+1)*8], c)
	}

	return L
}

// P - permutation of 64-bit block according to the table T.
// Permutation is necessary for the uniform distribution of bits in the block.
// The function takes a 64-bit block of bytes, and returns a new 64-bit block
// of bytes representing the permutation of the input block according to the
// table T.
func P(V []byte) []byte {
	T := [64]byte{
		0, 8, 16, 24, 32, 40, 48, 56, 1, 9, 17, 25, 33, 41, 49, 57,
		2, 10, 18, 26, 34, 42, 50, 58, 3, 11, 19, 27, 35, 43, 51, 59,
		4, 12, 20, 28, 36, 44, 52, 60, 5, 13, 21, 29, 37, 45, 53, 61,
		6, 14, 22, 30, 38, 46, 54, 62, 7, 15, 23, 31, 39, 47, 55, 63,
	}

	Result := make([]byte, 64)
	for i, t := range T {
		Result[i] = V[t]
	}
	return Result
}

// S is a substitution function that takes a 64-byte input and performs a
// permutation on it, returning a new 64-byte output. The permutation is defined
// by a 256-element permutation table, where each element is used to substitute
// the corresponding element in the input. The result is a new byte sequence
// that has the same length as the input.
func S(V []byte) []byte {
	Pi := [256]byte{
		252, 238, 221, 17, 207, 110, 49, 22, 251, 196, 250, 218, 35, 197, 4, 77,
		233, 119, 240, 219, 147, 46, 153, 186, 23, 54, 241, 187, 20, 205, 95, 193,
		249, 24, 101, 90, 226, 92, 239, 33, 129, 28, 60, 66, 139, 1, 142, 79,
		5, 132, 2, 174, 227, 106, 143, 160, 6, 11, 237, 152, 127, 212, 211, 31,
		235, 52, 44, 81, 234, 200, 72, 171, 242, 42, 104, 162, 253, 58, 206, 204,
		181, 112, 14, 86, 8, 12, 118, 18, 191, 114, 19, 71, 156, 183, 93, 135,
		21, 161, 150, 41, 16, 123, 154, 199, 243, 145, 120, 111, 157, 158, 178, 177,
		50, 117, 25, 61, 255, 53, 138, 126, 109, 84, 198, 128, 195, 189, 13, 87,
		223, 245, 36, 169, 62, 168, 67, 201, 215, 121, 214, 246, 124, 34, 185, 3,
		224, 15, 236, 222, 122, 148, 176, 188, 220, 232, 40, 80, 78, 51, 10, 74,
		167, 151, 96, 115, 30, 0, 98, 68, 26, 184, 56, 130, 100, 159, 38, 65,
		173, 69, 70, 146, 39, 94, 85, 47, 140, 163, 165, 125, 105, 213, 149, 59,
		7, 88, 179, 64, 134, 172, 29, 247, 48, 55, 107, 228, 136, 217, 231, 137,
		225, 27, 131, 73, 76, 63, 248, 254, 141, 83, 170, 144, 202, 216, 133, 97,
		32, 113, 103, 164, 45, 43, 9, 91, 203, 155, 37, 208, 190, 229, 108, 82,
		89, 166, 116, 210, 230, 244, 180, 192, 209, 102, 175, 194, 57, 75, 99, 182,
	}

	Result := make([]byte, 64)
	for i, b := range V {
		Result[i] = Pi[b]
	}
	return Result
}

// X performs a bitwise XOR operation on two 512-bit blocks represented as byte slices.
// The function takes two byte slices `k` and `a`, and returns a new byte slice
// containing the result of the XOR operation on each 64-bit segment of the input slices.

func X(k, a []byte) []byte {
	result := make([]byte, 64)
	for i := 0; i < 64; i += 8 {
		kVal := binary.LittleEndian.Uint64(k[i : i+8])
		aVal := binary.LittleEndian.Uint64(a[i : i+8])
		binary.LittleEndian.PutUint64(result[i:i+8], kVal^aVal)
	}
	return result
}

// K_array returns an array of 13 keys, each of length 64 bytes.
// The first key is the given key K_1, and the remaining keys are
// generated by applying the L, P, and S functions to the previous key
// and a set of fixed constants C[i]. The constants are hardcoded in the
// function. The resulting array of keys can be used to encrypt and decrypt
// data using the Streebog algorithm.

func K_array(K_1 []byte) [][]byte {
	C := [][]byte{
		mustHexDecode("b1085bda1ecadae9ebcb2f81c0657c1f2f6a76432e45d016714eb88d7585c4fc4b7ce09192676901a2422a08a460d31505767436cc744d23dd806559f2a64507"),
		mustHexDecode("6fa3b58aa99d2f1a4fe39d460f70b5d7f3feea720a232b9861d55e0f16b501319ab5176b12d699585cb561c2db0aa7ca55dda21bd7cbcd56e679047021b19bb7"),
		mustHexDecode("f574dcac2bce2fc70a39fc286a3d843506f15e5f529c1f8bf2ea7514b1297b7bd3e20fe490359eb1c1c93a376062db09c2b6f443867adb31991e96f50aba0ab2"),
		mustHexDecode("ef1fdfb3e81566d2f948e1a05d71e4dd488e857e335c3c7d9d721cad685e353fa9d72c82ed03d675d8b71333935203be3453eaa193e837f1220cbebc84e3d12e"),
		mustHexDecode("4bea6bacad4747999a3f410c6ca923637f151c1f1686104a359e35d7800fffbdbfcd1747253af5a3dfff00b723271a167a56a27ea9ea63f5601758fd7c6cfe57"),
		mustHexDecode("ae4faeae1d3ad3d96fa4c33b7a3039c02d66c4f95142a46c187f9ab49af08ec6cffaa6b71c9ab7b40af21f66c2bec6b6bf71c57236904f35fa68407a46647d6e"),
		mustHexDecode("f4c70e16eeaac5ec51ac86febf240954399ec6c7e6bf87c9d3473e33197a93c90992abc52d822c3706476983284a05043517454ca23c4af38886564d3a14d493"),
		mustHexDecode("9b1f5b424d93c9a703e7aa020c6e41414eb7f8719c36de1e89b4443b4ddbc49af4892bcb929b069069d18d2bd1a5c42f36acc2355951a8d9a47f0dd4bf02e71e"),
		mustHexDecode("378f5a541631229b944c9ad8ec165fde3a7d3a1b258942243cd955b7e00d0984800a440bdbb2ceb17b2b8a9aa6079c540e38dc92cb1f2a607261445183235adb"),
		mustHexDecode("abbedea680056f52382ae548b2e4f3f38941e71cff8a78db1fffe18a1b3361039fe76702af69334b7a1e6c303b7652f43698fad1153bb6c374b4c7fb98459ced"),
		mustHexDecode("7bcd9ed0efc889fb3002c6cd635afe94d8fa6bbbebab076120018021148466798a1d71efea48b9caefbacd1d7d476e98dea2594ac06fd85d6bcaa4cd81f32d1b"),
		mustHexDecode("378ee767f11631bad21380b00449b17acda43c32bcdf1d77f82012d430219f9b5d80ef9d1891cc86e71da4aa88e12852faf417d5d9b21b9948bc924af11bd720"),
	}

	Result := make([][]byte, 13)
	Result[0] = K_1

	for i := 1; i < 13; i++ {
		K_i := make([]byte, 64)
		for j := 0; j < 64; j += 8 {
			kVal := binary.LittleEndian.Uint64(Result[i-1][j : j+8])
			cVal := binary.LittleEndian.Uint64(C[i-1][j : j+8])
			binary.LittleEndian.PutUint64(K_i[j:j+8], kVal^cVal)
		}

		Result[i] = L(P(S(K_i)))
	}
	return Result
}

// mustHexDecode decodes a string containing a hexadecimal representation of a byte array,
// and panics if the decoding fails.
func mustHexDecode(s string) []byte {
	b, err := hex.DecodeString(s)
	if err != nil {
		panic(err)
	}
	return b
}

// E - encryption of a message block using a key expansion array
//
//	It takes an initial key K_1 and a message block m, and returns
//	the encrypted block. The encryption process involves a series of
//	linear, permutation, and substitution transformations using the
//	key expansion array.
func E(K_1 []byte, m []byte) []byte {
	K := K_array(K_1)
	var Result []byte = L(P(S(X(K[0], m))))
	for i := 1; i < 12; i++ {
		Result = L(P(S(X(K[i], Result))))
	}
	Result = X(K[12], Result)
	return Result
}

// GN - the hash function.
//
//	It takes an initial hash value h, a message block m, and a
//	nonce value N, and returns the updated hash value.
//	The function applies a series of linear, permutation, and
//	substitution transformations to the input values using the
//	key expansion array, and then XORs the result with the
//	initial hash value and the message block.
func GN(h, m, N []byte) []byte {
	temp := make([]byte, 64)

	// XOR(h, N)
	for i := 0; i < 64; i += 8 {
		hVal := binary.LittleEndian.Uint64(h[i : i+8])
		nVal := binary.LittleEndian.Uint64(N[i : i+8])
		binary.LittleEndian.PutUint64(temp[i:i+8], hVal^nVal)
	}

	K := L(P(S(temp)))
	E_m := E(K, m)

	Result := make([]byte, 64)

	// XOR(E_m, h) ⊕ m
	for i := 0; i < 64; i += 8 {
		eVal := binary.LittleEndian.Uint64(E_m[i : i+8])
		hVal := binary.LittleEndian.Uint64(h[i : i+8])
		mVal := binary.LittleEndian.Uint64(m[i : i+8])
		binary.LittleEndian.PutUint64(Result[i:i+8], eVal^hVal^mVal)
	}

	return Result
}

// Hash computes the hash of the given message using the GOST R 34.11-2012
// algorithm in the given mode (either 256 or 512). The mode parameter
// determines the length of the resulting hash. The function returns the
// resulting hash as a byte slice.
func Hash(Message []byte, mode int) []byte {
	var m []byte
	Mod := big.NewInt(0)
	Mod.Exp(big.NewInt(2), big.NewInt(512), nil)
	var IV []byte // Этап 1, инициализируем безобразие
	var h []byte
	IV = make([]byte, 64)
	var Sigma []byte = make([]byte, 64)
	var N []byte = make([]byte, 64)
	for i := 0; i < 64; i++ {
		IV[i] = 0x00
		Sigma[i] = 0x00
		N[i] = 0x00
	}

	if mode == 256 {
		for i := 0; i < 64; i++ {
			IV[i] = 0x01 // Заполняем каждый байт значением 0x01
		}
		// Теперь IV содержит 64 байта, каждый из которых равен 0x01
	}
	h = IV
	L_message := len(Message)
	for L_message >= 64 { // Начинается этап 2
		var M_help []byte = Message[:L_message-64]
		m = Message[L_message-64:]
		h = GN(h, m, N)
		arg_1 := new(big.Int)
		arg_2 := new(big.Int)
		var arg_3 *big.Int = big.NewInt(int64(512))
		arg_4 := new(big.Int)
		arg_1.SetBytes(Sigma)
		arg_2.SetBytes(m)
		arg_4.SetBytes(N)
		arg_3.Add(arg_4, arg_3)
		arg_3.Mod(arg_3, Mod)
		newN := make([]byte, len(N))
		arg_3.FillBytes(newN)
		N = newN
		arg_1.Add(arg_1, arg_2)
		arg_1.Mod(arg_1, Mod)
		newSigma := make([]byte, len(Sigma))
		arg_1.FillBytes(newSigma)
		Sigma = newSigma
		Message = M_help
		L_message = len(Message)
	}
	L_message = len(Message)

	// Начинается этап 3

	m = make([]byte, 64)
	index := 0
	for i := 0; i < 64-L_message-1; i++ {
		m[i] = 0x00
		index += 1
	}
	m[index] = 0x01
	index += 1
	for i := index; i < 64; i++ {
		m[i] = Message[L_message-64+i]
	}
	h = GN(h, m, N)
	arg_1 := new(big.Int)
	arg_2 := new(big.Int)
	var arg_3 *big.Int = big.NewInt(int64(L_message * 8))
	arg_4 := new(big.Int)
	arg_1 = arg_1.SetBytes(Sigma)
	arg_2 = arg_2.SetBytes(m)
	arg_4 = arg_4.SetBytes(N)
	arg_3.Add(arg_4, arg_3)
	arg_3.Mod(arg_3, Mod)
	newN := make([]byte, len(N))
	arg_3.FillBytes(newN)
	N = newN
	arg_1.Add(arg_1, arg_2)
	arg_1.Mod(arg_1, Mod)
	newSigma := make([]byte, len(Sigma))
	arg_1.FillBytes(newSigma)
	Sigma = newSigma
	O := make([]byte, 64)
	for i := 0; i < 64; i++ {
		O[i] = 0x00
	}
	h = GN(h, N, O)
	if mode == 256 {
		h = GN(h, Sigma, O)[0:32]
	} else if mode == 512 {
		h = GN(h, Sigma, O)
	}
	return h
}

// func main() {
// 	reader := bufio.NewReader(os.Stdin)

// 	fmt.Println("Вы хотите создать хеш вашего сообщения (1) или проверить целостность переданного сообщения (2)?")
// 	choice, _ := reader.ReadString('\n')
// 	choice = strings.TrimSpace(choice)

// 	switch choice {
// 	case "1":
// 		fmt.Println("Введите название файла с сообщением:")
// 		fileName, _ := reader.ReadString('\n')
// 		fileName = strings.TrimSpace(fileName)

// 		fmt.Println("Выберите режим работы (256 или 512):")
// 		modeStr, _ := reader.ReadString('\n')
// 		modeStr = strings.TrimSpace(modeStr)
// 		mode := 0
// 		if modeStr == "256" {
// 			mode = 256
// 		} else if modeStr == "512" {
// 			mode = 512
// 		} else {
// 			fmt.Println("Неверный режим работы. Допустимые значения: 256 или 512.")
// 			return
// 		}

// 		data, err := os.ReadFile(fileName)
// 		if err != nil {
// 			fmt.Println("Ошибка при чтении файла:", err)
// 			return
// 		}

// 		hash := Hash(data, mode)
// 		hashHex := hex.EncodeToString(hash)
// 		fmt.Println("Хеш сообщения:", hashHex)

// 		// Запись хеша в файл output_hash.txt
// 		err = os.WriteFile("output_hash.txt", []byte(hashHex), 0644)
// 		if err != nil {
// 			fmt.Println("Ошибка при записи хеша в файл:", err)
// 			return
// 		}
// 		fmt.Println("Хеш успешно записан в файл output_hash.txt")

// 	case "2":
// 		fmt.Println("Введите название файла с сообщением:")
// 		messageFile, _ := reader.ReadString('\n')
// 		messageFile = strings.TrimSpace(messageFile)

// 		fmt.Println("Введите название файла с хешем:")
// 		hashFile, _ := reader.ReadString('\n')
// 		hashFile = strings.TrimSpace(hashFile)

// 		fmt.Println("Выберите режим работы (256 или 512):")
// 		modeStr, _ := reader.ReadString('\n')
// 		modeStr = strings.TrimSpace(modeStr)
// 		mode := 0
// 		if modeStr == "256" {
// 			mode = 256
// 		} else if modeStr == "512" {
// 			mode = 512
// 		} else {
// 			fmt.Println("Неверный режим работы. Допустимые значения: 256 или 512.")
// 			return
// 		}

// 		messageData, err := os.ReadFile(messageFile)
// 		if err != nil {
// 			fmt.Println("Ошибка при чтении файла с сообщением:", err)
// 			return
// 		}

// 		hashData, err := os.ReadFile(hashFile)
// 		if err != nil {
// 			fmt.Println("Ошибка при чтении файла с хешем:", err)
// 			return
// 		}

// 		computedHash := Hash(messageData, mode)
// 		providedHash := strings.TrimSpace(string(hashData))

// 		if hex.EncodeToString(computedHash) == providedHash {
// 			fmt.Println("Целостность сообщения подтверждена.")
// 		} else {
// 			fmt.Println("Целостность сообщения нарушена.")
// 		}

// 	default:
// 		fmt.Println("Неверный выбор. Пожалуйста, выберите 1 или 2.")
// 	}
// }
