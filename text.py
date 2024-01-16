import numpy as np
def vectorise(mat):
    size = np.shape(mat)[0]
    vec = np.zeros(size**2,dtype=complex)
    for i in range(size):
        for j in range(size):
            vec[j+i*size] = mat[i,j]
    return vec
mat = np.array([[1+0j,2+0j],[3+0j,4+0j]])
print(mat)
vec = vectorise(mat)
print(vec)

def unvectorise(vec,size):
    mat = np.zeros((size,size))
    for k in range(size**2):
        i = k//size
        j = k%size
        mat[i,j] = vec[k]
    return mat

mat = unvectorise(vec,2)
print(mat)