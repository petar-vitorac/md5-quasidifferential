MODE = 1 # 0 = verify if same characteristic as Wang; 1 = get intermediate differences

# MD5 constants
K = [int(2**32*abs(sin(i+1))) for i in range(64)]
s = [7, 12, 17, 22]*4 + [5, 9, 14, 20]*4 + [4, 11, 16, 23]*4 + [6, 10, 15, 21]*4
a0 = int(0x67452301)
b0 = int(0xefcdab89)
c0 = int(0x98badcfe)
d0 = int(0x10325476)
BITS_PER_WORD = 32

# Differential from Wang's paper
dM = [0]*4 + [2**31] + [0]*6 + [2**15] + [0]*2 + [2**31] + [0]

#########################################
# UNCOMMENT THE MESSAGE YOU WISH TO USE #
#########################################

# Collision in Wang's paper
M0 = [0x2dd31d1, 0xc4eee6c5, 0x69a3d69, 0x5cf9af98, 0x87b5ca2f, 0xab7e4612, 0x3e580440, 0x897ffbb8, 0x634ad55, 0x2b3f409, 0x8388e483, 0x5a417125, 0xe8255108, 0x9fc9cdf7, 0xf2bd1dd9, 0x5b3c3780]

# Self computed collision 1
# M0 = [0xea2dc786, 0xc087389d, 0xc5f86464, 0xba921814, 0x14497da7, 0x5ecaf6d5, 0xe135826b, 0xd76678a6, 0x06352db1, 0x8473adf3, 0x8ba6b173, 0xd2b11966, 0xdff13cc3, 0x592f2f0d, 0x63768582, 0xf614e4a7, 0x670f844a, 0x76f77fd7, 0xfe55d9a5, 0xeb9f65a5, 0x25469175, 0x2cf4b284, 0x75ad0530, 0x8d9db84d, 0x6573382d, 0x86587727, 0x84198a42, 0x4c02bc6f, 0xc6224502, 0x5e7bcf94, 0x7e3e24cf, 0x233c973f]

# Self computed collision 2
# M0 = [0xe99819b4, 0x54920ae6, 0x98a1f04f, 0x521bd77a, 0x21f4fe1a, 0x6e79a1d2, 0x8a25903e, 0x4252c7b3, 0x0633ad41, 0x027377fb, 0x8a50e3b7, 0xe3515f6e, 0x946336d5, 0x804de62f, 0x68bbe190, 0xc5b8b922, 0x37ae5efa, 0xf79a3059, 0x1142eef9, 0xfb4e6222, 0x15f555d0, 0x1b8dd844, 0x665f5d36, 0xb0ad8ad5, 0x933af791, 0x67119543, 0x87176b92, 0xdceec084, 0x3c3dc8d1, 0xb55fc06f, 0x7eb42579, 0x5848e042]

# Self computed collision 3
# M0 = [0xb5d61dcc, 0x874f4bde, 0xfc3ef069, 0xd447be4a, 0xe6e2d622, 0x00976ecc, 0x5e1b87c4, 0x3203bdb7, 0x05336cd5, 0xc174940a, 0x7e7b005f, 0xd0d20b22, 0x34a50aef, 0x642e4ddc, 0xa9537b79, 0x30c5caf2, 0x3fdab570, 0x83b1b2ca, 0xa4d4f67b, 0xd38e1d0e, 0x1488d577, 0x1bf61840, 0x66bb6516, 0xc24144d7, 0x737b6ee9, 0x2e2652af, 0xb5744e3e, 0xc4429092, 0x584ab4c5, 0xddd3d054, 0x6defaa01, 0x9d2db194]

# Characteristic as in Wang's paper
char = [[]]*4 + [list(range(6,23)), [6,23,31], list(range(0,12))+list(range(23,32)), [0,23]+list(range(15,21)), [0,1,6,7,8,31], [12,13,31], [30,31],
    [7,8,31]+list(range(13,20)), [24,25,31], [31], [3,15,31], [29,31]] # which bits differ
char = [sum(map(lambda x: 2**x, el)) for el in char] # convert to integers

def leftrotate(x, n):
    return (x << n)|(x >> (BITS_PER_WORD-n))

# Executes MD5 on the given message and returns the intermediary results (those on the state registers if MODE=0, otherwise those between the state registers)
def md5(M):
    A = a0
    B = b0
    C = c0
    D = d0
    steps = list()
    for i in range(64):
        if i < 16:
            F = (B & C) | ((~B) & D)
            g = i
        elif i < 32:
            F = (D & B) | ((~D) & C)
            g = (5*i+1)%16
        elif i < 48:
            F = B ^^ C ^^ D
            g = (3*i+5)%16
        else:
            F = C ^^ (B | (~D))
            g = (7*i)%16
        
        if MODE == 1: # save intermediary results
            steps.append(int(F%(2**BITS_PER_WORD)))
            steps.append(int(F+A)%(2**BITS_PER_WORD))
            steps.append(int(F+A+M[g])%(2**BITS_PER_WORD))

        F = int((F + A + K[i] + M[g])%(2**BITS_PER_WORD))

        if MODE == 1: # save intermediary results
            steps.append(int(F))

        A = D
        D = C
        C = B
        B = int((B + leftrotate(F, s[i]))%(2**BITS_PER_WORD))

        if MODE == 0: # save results on the state registers
            steps.append(B)

    return steps

# M0' = M0 + dM (adding the differential to get the other message)
M0p = [M0i ^^ dMi for (M0i, dMi) in zip(M0, dM)]

H0 = md5(M0) # running MD5 on the message
H0p = md5(M0p) # running MD5 on the message + differential

observed_char = [H0i ^^ H0pi for (H0i, H0pi) in zip(H0, H0p)] # compute the differences by comparing the outputs

if MODE == 1:
    print("The intermediary differences are:")
    b1 = observed_char[::4]
    b2 = observed_char[1::4]
    b3 = observed_char[2::4]
    b4 = observed_char[3::4]
    print("b1 =", b1)
    print("b2 =", b2)
    print("b3 =", b3)
    print("b4 =", b4)
if MODE == 0:
    if not any([observed_char_i - char_i for (observed_char_i, char_i) in zip(observed_char, char)]):
        print("This message follows Wang's characteristic.")
    else:
        print("This message DOES NOT follow Wang's characteristic. The actual characteristic is:")
        print(observed_char)
