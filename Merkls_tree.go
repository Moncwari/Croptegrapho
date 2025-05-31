package main

import (
	"fmt"
	"os"
	"path/filepath"
)

func Merkl_tree_iteration(Values [][]byte) [][]byte {
	if len(Values)%2 == 1 {
		Values = append(Values, Values[len(Values)-1])
	}
	var Result [][]byte
	for i := 0; i < len(Values); i += 2 {
		combined := append(Values[i], Values[i+1]...)
		Result = append(Result, Hash(combined, 256))
		// Result = append(Result, Hash((Values[i]), 256))
		// Result = append(Result, Hash((Values[i+1]), 256))
	}
	return Result
}

func Merkl_tree(Values [][]byte) []byte {
	if len(Values) == 0 {
		return nil
	}
	for len(Values) > 1 {
		// fmt.Println(len(Values))
		Values = Merkl_tree_iteration(Values)
	}
	return Values[0]
}

func ReadRandomFiles(dir, pattern string) ([][]byte, error) {
	// Получаем список файлов
	files, err := filepath.Glob(filepath.Join(dir, pattern))
	if err != nil {
		return nil, fmt.Errorf("ошибка поиска файлов: %v", err)
	}

	if len(files) == 0 {
		return nil, fmt.Errorf("файлы не найдены (шаблон: %s)", pattern)
	}

	// Читаем файлы
	var result [][]byte
	for _, filename := range files {
		data, err := os.ReadFile(filename)
		if err != nil {
			return nil, fmt.Errorf("ошибка чтения файла %s: %v", filename, err)
		}
		result = append(result, data)
	}

	return result, nil
}

func getFirst5Bits(data []byte) byte {
	if len(data) == 0 {
		return 0
	}
	return (data[0] & 0xF8) >> 3 // Маска + сдвиг
}

func AppendToRandomFiles(data1, data2 []byte) error {
	for i := 1; i <= 5; i++ {
		filename := fmt.Sprintf("random_data_%d.bin", i)

		// Открываем файл для добавления (или создаем если не существует)
		file, err := os.OpenFile(filename, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
		if err != nil {
			return fmt.Errorf("ошибка открытия файла %s: %v", filename, err)
		}

		// Записываем первый []byte
		if _, err := file.Write(data1); err != nil {
			file.Close()
			return fmt.Errorf("ошибка записи в файл %s: %v", filename, err)
		}

		// Записываем второй []byte
		if _, err := file.Write(data2); err != nil {
			file.Close()
			return fmt.Errorf("ошибка записи в файл %s: %v", filename, err)
		}

		// Закрываем файл
		if err := file.Close(); err != nil {
			return fmt.Errorf("ошибка закрытия файла %s: %v", filename, err)
		}

		fmt.Printf("Данные дописаны в %s\n", filename)
	}
	return nil
}

// func main() {
// 	global := &Random_Numbers
// 	Files_byte, _ := ReadRandomFiles(".", "random_data_*.bin")
// 	// Merkl_Hash := Merkl_tree(Files_byte)
// 	//fmt.Printf("%x\n", Merkl_tree(Files_byte))
// 	//fmt.Println(BytesToHex(Merkl_tree(Files_byte)))
// 	//var All_signatures [][]byte
// 	for i := 1; i < 6; i++ {
// 		// Files_byte, _ := ReadRandomFiles(".", "random_data_*.bin")
// 		//fmt.Println("File №", i, ": ", BytesToHex(Hash(Files_byte[i], 256)))
// 		R, S := Schnorrs_signature(Files_byte[i-1])
// 		AppendToRandomFiles(R, S)
// 		// fmt.Println("Signature №", i, ": ", BytesToHex(R), BytesToHex(S))
// 		// All_signatures = append(All_signatures, signature)
// 	}
// 	Files_byte, _ = ReadRandomFiles(".", "random_data_*.bin")
// 	Merkl_Hash := Merkl_tree(Files_byte)

// 	First_part := []byte{0xFF, 0xFF, 0xFF, 0xFF} //make([]byte, 4)
// 	Second_part := Random_Number_Generator(*global, len(*global)+1)[len(*global)-1]
// 	Third_part := Merkl_Hash
// 	Fourth_part := []byte{
// 		0x01, // час
// 		0x1F, // день
// 		0x05, // месяц
// 		0x19, // год
// 	}
// 	nonce := make([]byte, 32)
// 	Heading := append(append(append(append(First_part, Second_part...), Third_part...), Fourth_part...), nonce...)
// 	Hash_find := Hash(Heading, 256)
// 	i := 0
// 	for getFirst5Bits(Hash_find) != 0 {
// 		increment(nonce)
// 		Heading = append(append(append(append(First_part, Second_part...), Third_part...), Fourth_part...), nonce...)
// 		Hash_find = Hash(Heading, 256)
// 		i++
// 	}
// 	fmt.Println(i)
// 	fmt.Println(Hash_find)
// 	// fmt.Println(BytesToHex(Merkl_tree(All_signatures)))
// }
