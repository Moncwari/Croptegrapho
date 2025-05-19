# Импорт твоего генератора и хэш-функции
from hashlib import new as streebog
import struct

# Используем генератор из предыдущего примера
class StreebogPRNG:
    def __init__(self, seed: bytes):
        self.state = stribog512(seed)

    def next_bytes(self, length: int) -> bytes:
        result = b''
        while len(result) < length:
            self.state = stribog512(self.state)
            result += self.state
        return result[:length]

    def next_int(self, bits: int = 256) -> int:
        byte_len = (bits + 7) // 8
        rand_bytes = self.next_bytes(byte_len)
        return int.from_bytes(rand_bytes, 'big') >> (byte_len * 8 - bits)

# Стрибог-512 (можно заменить на свою реализацию, если есть)
def stribog512(data: bytes) -> bytes:
    return streebog('streebog512', data).digest()

# Подпись Шнорра с key-prefixed, используя собственный ГПСЧ
def Schnorr(p, q, g, m, prng):
    r = prng.next_int(256) % q  # Нонс
    x = prng.next_int(256) % q  # Секретный ключ
    R = pow(g, r, p)
    P = pow(g, x, p)

    # Формируем сообщение для хэширования
    for_hash = R.to_bytes((R.bit_length() + 7) // 8, 'big') + \
               P.to_bytes((P.bit_length() + 7) // 8, 'big') + m
    e = int.from_bytes(stribog512(for_hash), 'big') % q
    s = (r + e * x) % q

    return p, q, g, P, R, s

# Проверка подписи Шнорра
def verif_Schnorr(p, q, g, P, R, s, m):
    for_hash = R.to_bytes((R.bit_length() + 7) // 8, 'big') + \
               P.to_bytes((P.bit_length() + 7) // 8, 'big') + m
    e = int.from_bytes(stribog512(for_hash), 'big') % q

    left = pow(g, s, p)
    right = (R * pow(P, e, p)) % p

    return left == right

# Пример использования
if __name__ == '__main__':
    p = int("EE8172AE8996608FB69359B89EB82A69854510E2977A4D63BC97322CE5DC3386EA0A12B343E9190F23177539845839786BB0C345D165976EF2195EC9B1C379E3", 16)
    q = int("98915E7EC8265EDFCDA31E88F24809DDB064BDC7285DD50D7289F0AC6F49DD2D", 16)
    g = int("9E96031500C8774A869582D4AFDE2127AFAD2538B4B6270A6F7C8837B50D50F206755984A49E509304D648BE2AB5AAB18EBE2CD46AC3D8495B142AA6CE23E21C", 16)

    seed = b'Student Lastname Firstname'  # замените на своё ФИО в байтах
    prng = StreebogPRNG(seed)

    text = "ГООООЛ"
    m = text.encode('utf-8')
    
    sign = Schnorr(p, q, g, m, prng)
    is_valid = verif_Schnorr(*sign[:5], sign[5], m)
    print("Подпись верна:", is_valid)
