from sage.crypto.sbox import SBox
S = SBox([7, 2, 4, 5, 1, 6, 3, 0])

def permute(x):
    y = vector(GF(2), 9)
    for i in range(9):
        s = i // 3 # S-box index
        j = i % 3  # Bit

        y[8 - (3 * j + ((s + 1) % 3))] = x[8 - i]
    return y

def sbox(x):
    y = vector(GF(2), 9)
    for i in range(0, 9, 3):
        y[i:i + 3] = S(x[i:i + 3])
    return y


def encrypt(x, rounds, roundkeys):
    y = copy(x)
    for i in range(rounds - 1):
        y = permute(sbox(y + roundkeys[i]))
    return sbox(y + roundkeys[rounds - 1]) + roundkeys[rounds]

roundkeys  = [vector(GF(2), [1, 0, 1, 1, 0, 0, 0, 1, 1]) for _ in range(3)] + [vector(GF(2), 9)]
# roundkeys  = [vector(GF(2), 9) for _ in range(4)]
plaintext  = vector(GF(2), [0, 0, 0, 0, 0, 0, 0, 0, 0])
ciphertext = encrypt(plaintext, 3, roundkeys)
print(ciphertext)

a = vector(GF(2), [0, 0, 0, 0, 0, 0, 0, 0, 1])
b = vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])

diff_trails = [
    [vector(GF(2), [0, 0, 1, 0, 0, 0, 0, 0, 0]), vector(GF(2), [0, 0, 0, 0, 0, 0, 0, 0, 1]), vector(GF(2), [0, 0, 0, 0, 0, 0, 0, 1, 0]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])],
    [vector(GF(2), [0, 0, 1, 0, 0, 0, 0, 0, 0]), vector(GF(2), [0, 0, 0, 0, 0, 0, 1, 0, 1]), vector(GF(2), [0, 1, 0, 0, 0, 0, 0, 1, 0]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])],
    [vector(GF(2), [0, 0, 1, 0, 0, 0, 0, 0, 0]), vector(GF(2), [0, 0, 0, 0, 0, 0, 0, 1, 1]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 1, 0]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])],
    [vector(GF(2), [0, 0, 1, 0, 0, 0, 0, 0, 0]), vector(GF(2), [0, 0, 0, 0, 0, 0, 1, 1, 1]), vector(GF(2), [0, 1, 0, 0, 1, 0, 0, 1, 0]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])]
]

total = 0

for differential in diff_trails:
    count = 0
    for plain in VectorSpace(GF(2), 9):
        ok = True
        for i in range(1, 4):
            ok = ok and encrypt(plain+a, i, roundkeys)==encrypt(plain, i, roundkeys)+differential[i]
        if ok:
            count += 1

    print(count/2**9)
    
    total += count
    
print("Total:")
print(total/2**9)

