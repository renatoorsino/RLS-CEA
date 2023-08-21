using LinearAlgebra
using Random

M = diagm([2, 1, 3, 2])
f = [1, 0, 1, 2]
g = 100
A = [g 0 0 0; 0 2*g 0 0; 1 -2 0 0]
b = [-g, -g, 0]

P = Dict()
S = Dict()
K = Dict()
H = Dict()
a = Dict()
e = Dict()

P[0] = inv(M)
a[0] = P[0]*f

s = shuffle(1:3)

for r in 0:2
    H[r+1] = A[s[r+1],:]'
    S[r+1] = H[r+1] * P[r] * H[r+1]' 
    if S[r+1] > 0 
        K[r+1] = P[r] * H[r+1]' / S[r+1]
        P[r+1] = P[r] - K[r+1] * S[r+1] * K[r+1]'
        e[r+1] = b[s[r+1]] - H[r+1] * a[r]
        a[r+1] = a[r] + K[r+1] * e[r+1]
    else
        P[r+1] = P[r]
        a[r+1] = a[r]
    end
end

print(s)
print("\n")
print(K)
print("\n")
print(e)
print("\n")
print(a[3])
print("\n")