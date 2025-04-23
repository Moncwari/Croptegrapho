import random
from gostcrypto import gosthash

def Schnorr(p, q, g, m):
    r = random.randint(0, q) # Нонс - определяем с ранее написанным ГСЧ
    x = random.randint(0, q) # Секретный ключ - определяем с ранее написанным ГСЧ
    R = pow(g, r, p)
    P = pow(g, x, p)

    for_hash = bytes.fromhex(hex(R)[2:]) + bytes.fromhex(hex(P)[2:]) + m
    hasher = gosthash.new('streebog512') # Хэш - определяем с ранее написанной функцией H
    hasher.update(for_hash)
    e = int.from_bytes(hasher.digest())
    s = (r + e*x) % q
    return p, q, g, P, R, s

def verif_Schnorr(p, q, g, P, R, s, m):
    for_hash = bytes.fromhex(hex(R)[2:]) + bytes.fromhex(hex(P)[2:]) + m
    hasher = gosthash.new('streebog512') # Хэш - определяем с ранее написанной функцией H
    hasher.update(for_hash)
    e = int.from_bytes(hasher.digest())

    #print(q)
    #print(pow(g, s, p))
    #print((R * pow(P, e, p)) % p)
    
    if (R * pow(P, e, p)) % p == pow(g, s, p):
        return True
    else:
        return False

p = int("EE8172AE8996608FB69359B89EB82A69854510E2977A4D63BC97322CE5DC3386EA0A12B343E9190F23177539845839786BB0C345D165976EF2195EC9B1C379E3", 16)
q = int("98915E7EC8265EDFCDA31E88F24809DDB064BDC7285DD50D7289F0AC6F49DD2D", 16)
g = int("9E96031500C8774A869582D4AFDE2127AFAD2538B4B6270A6F7C8837B50D50F206755984A49E509304D648BE2AB5AAB18EBE2CD46AC3D8495B142AA6CE23E21C", 16)

text = "ГООООЛ"
m = bytes(text, 'utf-8')
sign = Schnorr(p, q, g, m)
#print(sign)
print(verif_Schnorr(sign[0], sign[1], sign[2], sign[3], sign[4], sign[5], m))
