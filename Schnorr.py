from hash import Hash  # импорт твоей хэш-функции ГОСТ 34.11-2018
import time

# Класс генератора псевдослучайных чисел
class StreebogPRNG:
    def __init__(self, seed: bytes):
        self.state = Hash(seed, 512)  # начальное состояние

    def next_bytes(self, length: int) -> bytes:
        output = b''
        while len(output) < length:
            self.state = Hash(self.state, 512)
            output += self.state
        return output[:length]

    def next_int(self, bits: int = 256) -> int:
        byte_len = (bits + 7) // 8
        rnd = self.next_bytes(byte_len)
        return int.from_bytes(rnd, 'big') >> (byte_len * 8 - bits)

# Функция генерации подписи Шнорра
def Schnorr(p, q, g, m: bytes, prng: StreebogPRNG):
    r = prng.next_int(256) % q  # Нонс
    x = prng.next_int(256) % q  # Приватный ключ
    R = pow(g, r, p)
    P = pow(g, x, p)

    msg_hash_input = (
        R.to_bytes((R.bit_length() + 7) // 8, 'big') +
        P.to_bytes((P.bit_length() + 7) // 8, 'big') +
        m
    )
    e = int.from_bytes(Hash(msg_hash_input, 512), 'big') % q
    s = (r + e * x) % q
    return (p, q, g, P, R, s)

# Проверка подписи
def verif_Schnorr(p, q, g, P, R, s, m: bytes):
    msg_hash_input = (
        R.to_bytes((R.bit_length() + 7) // 8, 'big') +
        P.to_bytes((P.bit_length() + 7) // 8, 'big') +
        m
    )
    e = int.from_bytes(Hash(msg_hash_input, 512), 'big') % q
    left = pow(g, s, p)
    right = (R * pow(P, e, p)) % p
    return left == right

# Пример использования
if __name__ == "__main__":
    # Параметры ГОСТ Р 34.10–94, приложение А.3
    p = int("EE8172AE8996608FB69359B89EB82A69854510E2977A4D63BC97322CE5DC3386EA0A12B343E9190F23177539845839786BB0C345D165976EF2195EC9B1C379E3", 16)
    q = int("98915E7EC8265EDFCDA31E88F24809DDB064BDC7285DD50D7289F0AC6F49DD2D", 16)
    g = int("9E96031500C8774A869582D4AFDE2127AFAD2538B4B6270A6F7C8837B50D50F206755984A49E509304D648BE2AB5AAB18EBE2CD46AC3D8495B142AA6CE23E21C", 16)

    seed = b"Ivanov Ivan"  # замените на своё имя и фамилию
    prng = StreebogPRNG(seed)

    m = "ГООООЛ".encode("utf-8")
    sign = Schnorr(p, q, g, m, prng)
    valid = verif_Schnorr(*sign[:5], sign[5], m)

    print("Подпись корректна:", valid)
