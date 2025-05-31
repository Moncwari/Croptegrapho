package main

import (
	"encoding/binary"
	"fmt"
	"os"
)

func increment(bytes []byte) {
	for i := len(bytes) - 1; i >= 0; i-- {
		bytes[i]++
		if bytes[i] != 0 {
			break
		}
	}
}

//var I_start []byte
// var I_end []byte

func Random_Number_Generator(Random_Numbers [][]byte, I_end int) [][]byte {
	// Дополняем Message до 512 бит (64 байта)
	// paddedMessage := make([]byte, 64)
	// copy(paddedMessage, Message)

	// Нулевой цикл: H0 = Hash(Message)
	// H0 := Hash(paddedMessage, 256)
	var H0 []byte
	// var H_list [][]byte
	// H_list = append(H_list, H0)
	if len(Random_Numbers) == 0 {
		Message := []byte("Gorenkin Evgenii")
		paddedMessage := make([]byte, 64)
		copy(paddedMessage, Message)
		H0 := Hash(paddedMessage, 256)
		Random_Numbers = append(Random_Numbers, H0)
		// fmt.Println("check")
		if I_end == 1 {
			return Random_Numbers
		}
	}

	// else {
	// 	H0 = Random_Numbers[len(Random_Numbers)-1]
	// }

	I := make([]byte, 32)
	// binary.BigEndian.PutUint32(I, 1)
	// fmt.Println("Length: ", len(Random_Numbers))
	for i := len(Random_Numbers) - 1; i < I_end; i++ { // Генерируем 10 чисел (можно изменить)
		binary.BigEndian.PutUint32(I, uint32(i))
		// Формируем вход: H0 || I (8 + 4 = 12 байт? Нет, нужно 64!)
		// Надо дополнить I до 32 байт (256 бит), чтобы вместе с H0 было 512 бит.
		input := make([]byte, 64)
		copy(input[:32], H0)
		copy(input[32:], I)

		Hi := Hash(input, 256)
		// H_list = append(H_list, Hi)
		Random_Numbers = append(Random_Numbers, Hi)

		// Увеличиваем i на 1
		//increment(I)
		// print(I_start, "", I, "\n")
		// if bytes.Equal(I_start, I_end) {
		// 	break
		// }
	}

	return Random_Numbers //H_list
}

var (
	Random_Numbers [][]byte
)

func writeRandomDataToFiles() error {
	var Random_Numbers [][]byte
	totalBytesNeeded := 5 * 200                   // 5 файлов по 200 байт
	numbersNeeded := (totalBytesNeeded + 31) / 32 // Каждое число дает 32 байта

	// Генерируем все необходимые случайные числа за один вызов
	Random_Numbers = Random_Number_Generator(Random_Numbers, numbersNeeded)

	// Объединяем все числа в один большой байтовый массив
	var allData []byte
	for _, num := range Random_Numbers {
		allData = append(allData, num...)
	}

	// Записываем данные в файлы
	for i := 0; i < 5; i++ {
		start := i * 200
		end := start + 200
		if end > len(allData) {
			return fmt.Errorf("недостаточно сгенерированных данных")
		}

		filename := fmt.Sprintf("random_data_%d.bin", i+1)
		err := os.WriteFile(filename, allData[start:end], 0644)
		if err != nil {
			return fmt.Errorf("ошибка при записи файла %s: %v", filename, err)
		}
		fmt.Printf("Файл %s успешно создан (%d байт)\n", filename, 200)
	}

	return nil
}

// func main() {
// 	err := writeRandomDataToFiles()
// 	if err != nil {
// 		fmt.Println("Ошибка:", err)
// 	}
// }

// func main() {
// 	// name := []byte("Ivanov Ivan") // Пример входных данных
// 	//var Random_Numbers [][]byte
// 	// I := make([]byte, 32)

// 	for i := 0; i < 20; i++ {
// 		Random_Numbers = Random_Number_Generator(Random_Numbers, i)
// 		// fmt.Printf("%x\n", Random_Numbers[i-1])
// 		fmt.Println(i)
// 	}

// 	Random_Numbers = Random_Number_Generator(Random_Numbers, 20)
// 	fmt.Printf("%x\n", Random_Numbers[0])
// 	fmt.Println(len(Random_Numbers))
// }
