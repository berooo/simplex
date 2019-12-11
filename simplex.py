# -*- coding: utf-8 -*-

"""
@author: shibaorong
@contact: diamond_br@163.com
@Created on: 2019/12/8 12:12
"""
import numpy as np

def initialize(A,b,c):
    h,w=A.shape
    coefficient=np.zeros([h+1,w+1])
    bs=np.zeros([h+1,1])
    coefficient[0,:-1]=c
    coefficient[1:,:-1]=A
    coefficient[1:,-1]=b

    basis=np.asarray([x for x in range(w-h+1,w+1)])


    return coefficient,basis

def simplex(A,b,c):
    np.set_printoptions(suppress=True)
    """
    :param A:matrix
    :param b: vector
    :param c: vector
    :param mark: 表示c中有几个变量
    :return: optimal solution"""

    coeff,basis=initialize(A,b,c)

    h,w=coeff.shape

    negative_c_nums=-1
    last=-1

    while negative_c_nums!=0:
        negative_c_nums=0
        mark1=0
        mark_index1=0
        for i in range(w-1):
            if coeff[0,i]<mark1:
                negative_c_nums+=1
                mark1=coeff[0,i]
                mark_index1=i

        if negative_c_nums==0:
            break
        mark2=np.inf

        mark_index2=0
        for i in range(1,coeff.shape[0]):
            if coeff[i,mark_index1]>0:
                temp=coeff[i,w-1]/coeff[i,mark_index1]
                if temp==0 and temp==last:
                    continue
                if temp<mark2:
                    mark2=coeff[i,w-1]/coeff[i,mark_index1]
                    mark_index2=i
                else:
                    break
        last=temp
        #Gauss elimination
        coeff[mark_index2,:]/=coeff[mark_index2,mark_index1]
        for i in range(h):
            if coeff[i,mark_index1]!=0 and i!=mark_index2:
                coeff[i,:]-=coeff[mark_index2,:]*coeff[i,mark_index1]

        basis[mark_index2-1]=mark_index1+1

        np.set_printoptions(precision=2)

    z=-coeff[0,-1]
    X=np.zeros(coeff.shape[1]-1)
    for i in range(len(basis)):
        X[basis[i]-1]=coeff[i+1,-1]

    return z,X


if __name__=='__main__':

    A = np.array([[2, 3, 0, 1, 0, 0], [3, 5, -38, 0, 1, 0], [0, 0, 1, 0, 0, 1]])
    b = np.array([660, -600, 25])
    c = np.array([-60, -110, 510, 0, 0, 0])

    z,X=simplex(A,b,c)
    print(z)
    print(X)
