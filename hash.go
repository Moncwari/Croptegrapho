package main

import (
	"bufio"
	"errors"
	"fmt"
	"log"
	"math"
	"math/big"
	"os"
	"strconv"
	"strings"
)

// VecN converts an integer z into its binary representation
// as a string of n bits, padding with leading zeros if necessary.

func VecN(z int, n int) string {
	binary := fmt.Sprintf("%0*b", n, z) // Форматируем число в двоичное представление с ведущими нулями
	return binary
}

// BigVecN converts a big.Int z into its binary string representation
// of n bits. It pads the binary string with leading zeros if the
// length is less than n, or truncates the string from the left if
// it exceeds n bits.

func BigVecN(z *big.Int, n int) string {
	// Преобразуем число в двоичную строку
	binaryStr := z.Text(2)

	// Если длина строки меньше n, дополняем её ведущими нулями
	if len(binaryStr) < n {
		binaryStr = strings.Repeat("0", n-len(binaryStr)) + binaryStr
	}

	// Если длина строки больше n, обрезаем её до n символов
	if len(binaryStr) > n {
		binaryStr = binaryStr[len(binaryStr)-n:]
	}

	return binaryStr
}

// IntN converts a binary string of n bits into an integer.
// It panics if the length of the string is not equal to n, or if the string contains a character other than '0' or '1'.
// The output is the integer value of the binary string, with the least significant bit first (i.e., the rightmost bit is the least significant).
func IntN(binaryStr string, n int) int {
	// Проверяем, что длина строки соответствует n
	if len(binaryStr) != n {
		panic("Длина двоичной строки не соответствует n")
	}

	// Инициализируем результат
	var z uint

	// Проходим по каждому биту строки
	for i := 0; i < n; i++ {
		// Получаем текущий бит (символ '0' или '1')
		bit := binaryStr[n-1-i] // Переворачиваем порядок битов

		// Преобразуем символ в число (0 или 1)
		if bit == '1' {
			z += uint(math.Pow(2, float64(i)))
		} else if bit != '0' {
			fmt.Println(bit, string(binaryStr[n-1-i]))
			panic("Недопустимый символ в двоичной строке")
		}
	}

	return int(z)
}

// BigIntN converts a binary string of n bits into a big.Int.
// It panics if the length of the string is not equal to n, or if the string contains a character other than '0' or '1'.
// The output is the big.Int value of the binary string, with the least significant bit first (i.e., the rightmost bit is the least significant).
func BigIntN(binaryStr string, n int) *big.Int {
	// Проверяем, что длина строки соответствует n
	if len(binaryStr) != n {
		panic("Длина двоичной строки не соответствует n")
	}

	// Инициализируем результат как big.Int
	result := new(big.Int)

	// Проходим по каждому биту строки
	for i := 0; i < n; i++ {
		// Получаем текущий бит (символ '0' или '1')
		bit := binaryStr[n-1-i] // Переворачиваем порядок битов

		// Преобразуем символ в число (0 или 1)
		if bit == '1' {
			// Добавляем 2^i к результату
			result.SetBit(result, i, 1)
		} else if bit != '0' {
			panic("Недопустимый символ в двоичной строке")
		}
	}

	return result
}

// L performs a linear transformation on a string V of 512 bits.
// The result is a string of up to 512 bits.
// The transformation is performed by multiplying each block of 64 bits of V
// by a 4x64 matrix, and then XORing the result with the original block.
func L(V string) string { // Линейное преобразование
	var L string
	var Matrix [16][4]string = [16][4]string{
		{"8e20faa72ba0b470", "47107ddd9b505a38", "ad08b0e0c3282d1c", "d8045870ef14980e"},
		{"6c022c38f90a4c07", "3601161cf205268d", "1b8e0b0e798c13c8", "83478b07b2468764"},
		{"a011d380818e8f40", "5086e740ce47c920", "2843fd2067adea10", "14aff010bdd87508"},
		{"0ad97808d06cb404", "05e23c0468365a02", "8c711e02341b2d01", "46b60f011a83988e"},
		{"90dab52a387ae76f", "486dd4151c3dfdb9", "24b86a840e90f0d2", "125c354207487869"},
		{"092e94218d243cba", "8a174a9ec8121e5d", "4585254f64090fa0", "accc9ca9328a8950"},
		{"9d4df05d5f661451", "c0a878a0a1330aa6", "60543c50de970553", "302a1e286fc58ca7"},
		{"18150f14b9ec46dd", "0c84890ad27623e0", "0642ca05693b9f70", "0321658cba93c138"},
		{"86275df09ce8aaa8", "439da0784e745554", "afc0503c273aa42a", "d960281e9d1d5215"},
		{"e230140fc0802984", "71180a8960409a42", "b60c05ca30204d21", "5b068c651810a89e"},
		{"456c34887a3805b9", "ac361a443d1c8cd2", "561b0d22900e4669", "2b838811480723ba"},
		{"9bcf4486248d9f5d", "c3e9224312c8c1a0", "effa11af0964ee50", "f97d86d98a327728"},
		{"e4fa2054a80b329c", "727d102a548b194e", "39b008152acb8227", "9258048415eb419d"},
		{"492c024284fbaec0", "aa16012142f35760", "550b8e9e21f7a530", "a48b474f9ef5dc18"},
		{"70a6a56e2440598e", "3853dc371220a247", "1ca76e95091051ad", "0edd37c48a08a6d8"},
		{"07e095624504536c", "8d70c431ac02a736", "c83862965601dd1b", "641c314b2b8ee083"}}
	for l := 0; l < 8; l++ {
		var c string = strings.Repeat("0", 64)
		b := V[l*64 : (l+1)*64]
		for i := 0; i < 64; i++ {
			var c_help string
			b_help := string(b[i])
			var vec_help_1 string
			vec_help_1 = Matrix[i/4][i%4]
			vec_help_1 = stringToBinary(vec_help_1)
			if b_help == "0" {
				c_help += strings.Repeat("0", 64)
			} else if b_help == "1" {
				c_help += vec_help_1
			}
			a, _ := strconv.ParseUint(c, 2, 64)
			b, _ := strconv.ParseUint(c_help, 2, 64)
			c = strconv.FormatUint(a^b, 2)
			if len(c) < 64 {
				c = strings.Repeat("0", 64-len(c)) + c
			}
		}
		if c != "0" {
			L += c
		}
	}
	return L
}

// P - permutation of 64-bit block according to the table T.
// Permutation is necessary for the uniform distribution of bits in the block.
func P(V string) string { //Перестановка
	var Result string
	var T [64]int = [64]int{0, 8, 16, 24, 32, 40, 48, 56, 1, 9, 17, 25, 33, 41, 49, 57, 2, 10, 18, 26, 34, 42, 50, 58,
		3, 11, 19, 27, 35, 43, 51, 59, 4, 12, 20, 28, 36, 44, 52, 60, 5, 13, 21, 29, 37, 45, 53, 61, 6, 14,
		22, 30, 38, 46, 54, 62, 7, 15, 23, 31, 39, 47, 55, 63}
	for i := 0; i < 64; i++ {
		Result += V[T[i]*8 : (T[i]+1)*8]
	}
	return Result
}

// S - nonlinear transformation of 64-bit block according to the table Pi
//
//	It is used in the key generation algorithm.
func S(V string) string {
	var Result string
	var Pi_help [256]int = [256]int{252, 238, 221, 17, 207, 110, 49, 22, 251, 196, 250, 218, 35, 197, 4, 77, 233, 119, 240,
		219, 147, 46, 153, 186, 23, 54, 241, 187, 20, 205, 95, 193, 249, 24, 101, 90, 226, 92, 239,
		33, 129, 28, 60, 66, 139, 1, 142, 79, 5, 132, 2, 174, 227, 106, 143, 160, 6, 11, 237, 152, 127,
		212, 211, 31, 235, 52, 44, 81, 234, 200, 72, 171, 242, 42, 104, 162, 253, 58, 206, 204, 181,
		112, 14, 86, 8, 12, 118, 18, 191, 114, 19, 71, 156, 183, 93, 135, 21, 161, 150, 41, 16, 123,
		154, 199, 243, 145, 120, 111, 157, 158, 178, 177, 50, 117, 25, 61, 255, 53, 138, 126, 109,
		84, 198, 128, 195, 189, 13, 87, 223, 245, 36, 169, 62, 168, 67, 201, 215, 121, 214, 246, 124,
		34, 185, 3, 224, 15, 236, 222, 122, 148, 176, 188, 220, 232, 40, 80, 78, 51, 10, 74, 167, 151,
		96, 115, 30, 0, 98, 68, 26, 184, 56, 130, 100, 159, 38, 65, 173, 69, 70, 146, 39, 94, 85, 47,
		140, 163, 165, 125, 105, 213, 149, 59, 7, 88, 179, 64, 134, 172, 29, 247, 48, 55, 107, 228,
		136, 217, 231, 137, 225, 27, 131, 73, 76, 63, 248, 254, 141, 83, 170, 144, 202, 216, 133, 97,
		32, 113, 103, 164, 45, 43, 9, 91, 203, 155, 37, 208, 190, 229, 108, 82, 89, 166, 116, 210,
		230, 244, 180, 192, 209, 102, 175, 194, 57, 75, 99, 182}

	for i := 0; i < 64; i++ {
		Result += VecN(Pi_help[IntN(V[i*8:(i+1)*8], 8)], 8)
	}
	return Result
}

// X - bit-by-bit XOR operation of two 512-bit blocks
//
//	It is used in the key generation algorithm and in the hash function.
func X(k string, a string) string {

	K_help := new(big.Int)
	A_help := new(big.Int)
	_, err := K_help.SetString(k, 2)
	_, err = A_help.SetString(a, 2)
	if err {
	}
	arg := new(big.Int)
	arg.Xor(K_help, A_help)
	answer := arg.Text(2)
	if len(answer) < 512 {
		answer = strings.Repeat("0", 512-len(answer)) + arg.Text(2)
	}
	return answer
}

// K_array generates an array of 13 strings used for key expansion.
// It takes an initial key K_1 and performs transformations using a predefined
// set of constants C. The resulting array is used in cryptographic operations.
// Each subsequent key is derived by XORing the previous key with a constant,
// followed by a series of linear, permutation, and substitution transformations.

func K_array(K_1 string) [13]string {
	var C [12]string
	var C_help [12]string = [12]string{"b1085bda1ecadae9ebcb2f81c0657c1f2f6a76432e45d016714eb88d7585c4fc4b7ce09192676901a2422a08a460d31505767436cc744d23dd806559f2a64507", "6fa3b58aa99d2f1a4fe39d460f70b5d7f3feea720a232b9861d55e0f16b501319ab5176b12d699585cb561c2db0aa7ca55dda21bd7cbcd56e679047021b19bb7", " f574dcac2bce2fc70a39fc286a3d843506f15e5f529c1f8bf2ea7514b1297b7bd3e20fe490359eb1c1c93a376062db09c2b6f443867adb31991e96f50aba0ab2", "ef1fdfb3e81566d2f948e1a05d71e4dd488e857e335c3c7d9d721cad685e353fa9d72c82ed03d675d8b71333935203be3453eaa193e837f1220cbebc84e3d12e", "4bea6bacad4747999a3f410c6ca923637f151c1f1686104a359e35d7800fffbdbfcd1747253af5a3dfff00b723271a167a56a27ea9ea63f5601758fd7c6cfe57", "ae4faeae1d3ad3d96fa4c33b7a3039c02d66c4f95142a46c187f9ab49af08ec6cffaa6b71c9ab7b40af21f66c2bec6b6bf71c57236904f35fa68407a46647d6e", "f4c70e16eeaac5ec51ac86febf240954399ec6c7e6bf87c9d3473e33197a93c90992abc52d822c3706476983284a05043517454ca23c4af38886564d3a14d493", "9b1f5b424d93c9a703e7aa020c6e41414eb7f8719c36de1e89b4443b4ddbc49af4892bcb929b069069d18d2bd1a5c42f36acc2355951a8d9a47f0dd4bf02e71e", "378f5a541631229b944c9ad8ec165fde3a7d3a1b258942243cd955b7e00d0984800a440bdbb2ceb17b2b8a9aa6079c540e38dc92cb1f2a607261445183235adb", "abbedea680056f52382ae548b2e4f3f38941e71cff8a78db1fffe18a1b3361039fe76702af69334b7a1e6c303b7652f43698fad1153bb6c374b4c7fb98459ced", "7bcd9ed0efc889fb3002c6cd635afe94d8fa6bbbebab076120018021148466798a1d71efea48b9caefbacd1d7d476e98dea2594ac06fd85d6bcaa4cd81f32d1b", "378ee767f11631bad21380b00449b17acda43c32bcdf1d77f82012d430219f9b5d80ef9d1891cc86e71da4aa88e12852faf417d5d9b21b9948bc924af11bd720"}
	var Result [13]string
	for i := 0; i < 12; i++ {
		C[i] = stringToBinary(C_help[i])

	}
	Result[0] = K_1
	for i := 1; i < 13; i++ {
		K_help := new(big.Int)
		C_help := new(big.Int)
		_, err := K_help.SetString(Result[i-1], 2)
		_, err = C_help.SetString(C[i-1], 2)
		if err {
		}
		arg := new(big.Int)
		arg.Xor(K_help, C_help)
		answer := arg.Text(2)
		if len(answer) < 512 {
			answer = strings.Repeat("0", 512-len(answer)) + arg.Text(2)
		}
		Result[i] = L(P(S(answer)))
	}
	return Result
}

// E - encryption of a message block using a key expansion array
//
//	It takes an initial key K_1 and a message block m, and returns
//	the encrypted block. The encryption process involves a series of
//	linear, permutation, and substitution transformations using the
//	key expansion array.
func E(K_1 string, m string) string {
	K := K_array(K_1)
	var Result string = L(P(S(X(K[0], m))))
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
func GN(h string, m string, N string) string {
	var Result string
	h_help := new(big.Int)
	N_help := new(big.Int)
	m_help := new(big.Int)
	_, err := h_help.SetString(h, 2)
	_, err = N_help.SetString(N, 2)
	_, err = m_help.SetString(m, 2)
	if err {
	}
	arg := new(big.Int)
	arg.Xor(h_help, N_help)
	answer := arg.Text(2)
	for len(answer) < 512 {
		answer = "0" + answer
	}
	Result = E(L(P(S(answer))), m)
	Result_help := new(big.Int)
	_, err = Result_help.SetString(Result, 2)
	Result_help.Xor(Result_help, h_help)
	Result = Result_help.Xor(Result_help, m_help).Text(2)
	if len(Result) < 512 {
		Result = strings.Repeat("0", 512-len(Result)) + Result
	}
	return Result
}

// binaryToHex converts a binary string of any length into a hexadecimal string.
// The function takes care to ensure that the binary string length is a multiple of 4,
// by adding leading zeros if necessary. It then processes the string four bits at
// a time, converting each 4-bit chunk into a hexadecimal digit, and adding it
// to the result string. The function ignores any errors that may occur during the
// conversion process.
func binaryToHex(binaryString string) string {
	// Убедимся, что длина строки кратна 4, добавив ведущие нули при необходимости
	for len(binaryString)%4 != 0 {
		binaryString = "0" + binaryString
	}

	// Результирующая шестнадцатеричная строка
	hexString := ""

	// Проходим по строке с шагом 4
	for i := 0; i < len(binaryString); i += 4 {
		// Берем 4 бита
		chunk := binaryString[i : i+4]

		// Преобразуем 4 бита в число (0-15)
		value, err := strconv.ParseUint(chunk, 2, 64)
		if err != nil {
			// В случае ошибки просто пропускаем этот фрагмент
			continue
		}

		// Преобразуем число в шестнадцатеричную цифру
		hexDigit := fmt.Sprintf("%X", value)

		// Добавляем цифру к результату
		hexString += hexDigit
	}

	return hexString
}

// stringToBinary преобразует строку в 16-ричном формате в бинарную строку.
// Функция игнорирует любые ошибки, которые могут возникнуть при преобразовании,
// просто пропуская невалидные символы.
//
// hexString - строка в 16-ричном формате, которую нужно преобразовать.
// Возвращает строку, состоящую из 0 и 1, представляющую бинарное представление
// входной строки.
func stringToBinary(hexString string) string {
	// Убираем возможные пробелы или префиксы
	hexString = strings.ReplaceAll(hexString, " ", "")
	hexString = strings.TrimPrefix(hexString, "0x")

	// Результирующая бинарная строка
	binaryString := ""

	// Проходим по каждому символу в строке
	for _, char := range hexString {
		// Преобразуем символ в число (0-15)
		var value int
		if char >= '0' && char <= '9' {
			value = int(char - '0')
		} else if char >= 'A' && char <= 'F' {
			value = int(char-'A') + 10
		} else if char >= 'a' && char <= 'f' {
			value = int(char-'a') + 10
		} else {
			// Если символ невалидный, пропускаем его
			continue
		}

		// Преобразуем число в 4-битное бинарное представление
		binaryString += fmt.Sprintf("%04b", value)
	}

	return binaryString
}

// Hash calculates the hash of a given Message string using the GOST R 34.11-2012 (Streebog) hash function.
// The function takes a Message string and a mode (either 256 or 512) and returns the calculated hash as a string.
// The hash calculation involves padding the message to a multiple of 512 bits, then applying a series of linear,
// permutation, and substitution transformations to the padded message and the hash value, using the key expansion array.
// The resulting hash value is then returned as a string.
func Hash(Message string, mode int) string {
	var m string
	Mod := big.NewInt(0)
	Mod.Exp(big.NewInt(2), big.NewInt(512), nil)
	M := stringToBinary(Message)
	var IV string // Этап 1, инициализируем безобразие
	var h string
	for i := 0; i < 512; i++ {
		IV += "0"
	}
	var Sigma string = IV
	var N string = IV
	if mode == 256 {
		IV = ""
		for i := 0; i < 64; i++ {
			IV += "00000001"
		}
	}
	h = IV

	for len(M) >= 512 { // Начинается этап 2
		var M_help string = M[:len(M)-512]
		m = M[len(M)-512:]
		h = GN(h, m, N)
		arg_1 := BigIntN(Sigma, 512)
		arg_2 := BigIntN(m, 512)
		arg_3 := big.NewInt(int64(512))
		arg_4 := BigIntN(N, 512)
		arg_3.Add(arg_4, arg_3)
		arg_3.Mod(arg_3, Mod)
		N = BigVecN(arg_3, 512)
		arg_1.Add(arg_1, arg_2)
		arg_1.Mod(arg_1, Mod)
		Sigma = BigVecN(arg_1, 512)
		M = M_help
	}
	m = ""
	for i := 0; i < 511-len(M); i++ { // Начинается этап 3
		m += "0"
	}
	m += "1"
	m += M
	h = GN(h, m, N)
	arg_1 := BigIntN(Sigma, 512)
	arg_2 := BigIntN(m, 512)
	arg_3 := big.NewInt(int64(len(M)))
	arg_4 := BigIntN(N, 512)
	arg_3.Add(arg_4, arg_3)
	arg_3.Mod(arg_3, Mod)
	N = BigVecN(arg_3, 512)
	arg_1.Add(arg_1, arg_2)
	arg_1.Mod(arg_1, Mod)
	Sigma = BigVecN(arg_1, 512)
	h = GN(h, N, "0")
	if mode == 256 {
		h = GN(h, Sigma, "0")[0:256]
	} else if mode == 512 {
		h = GN(h, Sigma, "0")
	}
	return h
}

func IncrementalHash(mode int) (func(string), func() string) {
	var IV string
	var h string
	var Sigma string
	var N string
	Mod := big.NewInt(0)
	Mod.Exp(big.NewInt(2), big.NewInt(512), nil)

	// Инициализация состояния
	for i := 0; i < 512; i++ {
		IV += "0"
	}
	Sigma = IV
	N = IV
	if mode == 256 {
		IV = ""
		for i := 0; i < 64; i++ {
			IV += "00000001"
		}
	}
	h = IV

	// Функция для обновления состояния хэширования
	update := func(m string) {
		// Дополняем блок до 512 бит, если его длина меньше
		if len(m) < 512 {
			m = m + "1" + strings.Repeat("0", 511-len(m))
		}
		arg_1 := BigIntN(Sigma, 512)
		arg_2 := BigIntN(m, 512)
		arg_3 := big.NewInt(int64(len(m)))
		arg_4 := BigIntN(N, 512)
		arg_3.Add(arg_4, arg_3)
		arg_3.Mod(arg_3, Mod)
		N = BigVecN(arg_3, 512)
		arg_1.Add(arg_1, arg_2)
		arg_1.Mod(arg_1, Mod)
		Sigma = BigVecN(arg_1, 512)
		h = GN(h, m, N)
	}

	// Функция для финализации хэширования и получения итогового хэша
	finalize := func() string {
		h = GN(h, N, "0")
		if mode == 256 {
			h = GN(h, Sigma, "0")[0:256]
		} else if mode == 512 {
			h = GN(h, Sigma, "0")
		}
		return h
	}

	return update, finalize
}

// HashFile читает файл по блокам и вычисляет хэш
func HashFile(inputFile, outputFile string, mode int) error {
	// Открываем файл для чтения
	file, err := os.Open(inputFile)
	if err != nil {
		return err
	}
	defer file.Close()

	// Инициализируем инкрементальное хэширование
	update, finalize := IncrementalHash(mode)

	// Читаем файл по блокам
	buffer := make([]byte, 64) // 512 бит = 64 байта
	fileSize := 0
	for {
		n, err := file.Read(buffer)
		if err != nil && n == 0 {
			break
		}
		fileSize += n
		if n > 0 {
			// Преобразуем байты в бинарную строку
			binaryString := ""
			for _, b := range buffer[:n] {
				binaryString += fmt.Sprintf("%08b", b)
			}
			update(binaryString)
		}
	}

	// Дополняем последний блок, если его размер меньше 512 бит
	if fileSize%64 != 0 {
		lastBlock := make([]byte, 64)
		copy(lastBlock, buffer[:fileSize%64])
		lastBlock[fileSize%64] = 0x80 // Добавляем бит "1"
		for i := fileSize%64 + 1; i < 64; i++ {
			lastBlock[i] = 0 // Дополняем нулями
		}
		binaryString := ""
		for _, b := range lastBlock {
			binaryString += fmt.Sprintf("%08b", b)
		}
		update(binaryString)
	}

	// Добавляем блок с длиной сообщения
	lengthBlock := make([]byte, 64)
	length := uint64(fileSize * 8) // Длина в битах
	for i := 0; i < 8; i++ {
		lengthBlock[63-i] = byte(length >> (8 * i))
	}
	binaryString := ""
	for _, b := range lengthBlock {
		binaryString += fmt.Sprintf("%08b", b)
	}
	update(binaryString)

	// Финализируем хэширование
	hash := finalize()

	// Записываем результат в файл
	output, err := os.Create(outputFile)
	if err != nil {
		return err
	}
	defer output.Close()

	_, err = output.WriteString(binaryToHex(hash))
	if err != nil {
		return err
	}

	return nil
}

// ReadFileAsString reads the contents of a file and returns it as a string.
// If an error occurs during the read, an empty string and the error are returned.
// The function is intended for reading small files, such as key files or other
// configuration files. For larger files, use a streaming approach to avoid
// loading the entire file into memory.
func ReadFileAsString(filename string) (string, error) {
	data, err := os.ReadFile(filename)
	if err != nil {
		return "", err
	}
	return string(data), nil
}

// WriteStringToFile writes the given string to a file. If the file does not exist,
// it is created. If the file does exist, its contents are overwritten. The
// function returns an error if the file cannot be written to, e.g. if the
// directory does not exist, or if the user does not have the necessary
// permissions.
func WriteStringToFile(filename, data string) error {
	return os.WriteFile(filename, []byte(data), 0644)
}

// readChoice2 asks the user to enter a choice (1 or 2) and returns the input as an integer.
// The function returns an error if the input is not a valid integer or is not 1 or 2.
// The function is used to ask the user whether they want to check the integrity of a file
// or calculate the hash of a file.
func readChoice2() (int, error) {
	reader := bufio.NewReader(os.Stdin)

	fmt.Print("Введите 1 или 2: ")
	input, err := reader.ReadString('\n')
	if err != nil {
		return 0, fmt.Errorf("ошибка чтения ввода: %w", err)
	}

	input = strings.TrimSpace(input)

	number, err := strconv.Atoi(input)
	if err != nil {
		return 0, fmt.Errorf("некорректный формат числа: %w", err)
	}

	if number != 1 && number != 2 {
		return 0, errors.New("введено некорректное число. Ожидается 1 или 2")
	}

	return number, nil
}

func main() {
	//M := "fbe2e5f0eee3c820fbeafaebef20fffbf0e1e0f0f520e0ed20e8ece0ebe5f0f2f120fff0eeec20f120faf2fee5e2202ce8f6f3ede220e8e6eee1e8f0f2d1202ce8f0f2e5e220e5d1" //"323130393837363534333231303938373635343332313039383736353433323130393837363534333231303938373635343332313039383736353433323130" //"fbe2e5f0eee3c820fbeafaebef20fffbf0e1e0f0f520e0ed20e8ece0ebe5f0f2f120fff0eeec20f120faf2fee5e2202ce8f6f3ede220e8e6eee1e8f0f2d1202ce8f0f2e5e220e5d1" //
	//fmt.Println(binaryToHex(Hash(M, 512)))
	fmt.Println("Вы хотите проверить целостность файла или вычислить хэш? (1 - проверка, 2 - вычисление)")
	choice_1, _ := readChoice2()
	if choice_1 == 1 {
	}
	fmt.Println("Введите путь к файлу с хэш-кодом (например, hash.txt.txt):")
	var hash string
	hash, _ = ReadFileAsString("hash.txt")
	var choice_3 string
	fmt.Scanln(&choice_3)
	fmt.Println("Введите путь к файлу и режим работы через пробел (например, input.txt 512):")
	var choice_2 string
	fmt.Scanln(&choice_2)
	filename := "input.txt" // Укажите путь к файлу
	content, err := ReadFileAsString(filename)
	//fmt.Println(content)
	if err != nil {
		log.Fatalf("Ошибка при чтении файла: %v", err)
	}
	answer := binaryToHex(Hash(content, 256))
	if answer == hash {
		fmt.Println("Файл целостен")
	} else {
		fmt.Println("Файл поврежден")
	}
	//fmt.Println(answer, len(Hash(content, 512)))
	WriteStringToFile("output_hash.txt", answer)
	//fmt.Println("Хэш успешно вычислен и сохранен в output_hash.txt")

	// inputFile := "input.txt"
	// outputFile := "output_hash.txt"
	// mode := 512 // или 256

	// err := HashFile(inputFile, outputFile, mode)
	// if err != nil {
	// 	fmt.Println("Ошибка:", err)
	// } else {
	// 	fmt.Println("Хэш успешно вычислен и сохранен в", outputFile)
	// }
}
