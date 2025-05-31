package main

import "fmt"

func main() {
	err := writeRandomDataToFiles()
	if err != nil {
		fmt.Println("Ошибка:", err)
	}

	global := &Random_Numbers
	Files_byte, _ := ReadRandomFiles(".", "random_data_*.bin")
	for i := 1; i < 6; i++ {
		R, S := Schnorrs_signature(Files_byte[i-1])
		AppendToRandomFiles(R, S)
	}
	Files_byte, _ = ReadRandomFiles(".", "random_data_*.bin")
	Merkl_Hash := Merkl_tree(Files_byte)

	First_part := []byte{0xFF, 0xFF, 0xFF, 0xFF} //make([]byte, 4)
	Second_part := Random_Number_Generator(*global, len(*global)+1)[len(*global)-1]
	Third_part := Merkl_Hash
	Fourth_part := []byte{
		0x01, // час
		0x1F, // день
		0x05, // месяц
		0x19, // год
	}
	nonce := make([]byte, 32)
	Heading := append(append(append(append(First_part, Second_part...), Third_part...), Fourth_part...), nonce...)
	Hash_find := Hash(Heading, 256)
	i := 0
	for getFirst5Bits(Hash_find) != 0 {
		increment(nonce)
		Heading = append(append(append(append(First_part, Second_part...), Third_part...), Fourth_part...), nonce...)
		Hash_find = Hash(Heading, 256)
		i++
	}
	fmt.Println(i)
	fmt.Println(Hash_find)
}
