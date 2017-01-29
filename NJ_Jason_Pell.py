#!/usr/bin/python

'''
*******************************************************************************
tree.py

By: Jason Pell

An implementation of the neighbor-joining algorithm.

References:
N Saitou, M Nei. The neighbor joining method: a new method for reconstructing
phylogenetic trees. Molecular Biology and Evolution 4 (1987): no. 4, 406-425.
*******************************************************************************
'''

import numpy

def calcDist(distMat, i, j):
   ''' 
   Returns the distance between two taxa.

   Receive: Distance matrix and the two taxa.
   Return: Distance between taxa.
   '''

   if i < j:
      i, j = j, i
   return distMat[i][j]

def calcDistSum(distMat, i):
   '''
   Calculates the sum of distances for a taxa.

   Receive: Distance matrix and a taxa.
   Return: The sum of distances.
   '''

   sum = 0

   for k in range(len(distMat)):
      sum += distMat[i][k]

   for k in range(len(distMat)):
      sum += distMat[k][i]

   return sum

def calcQ(distMat, i, j):
   '''
   Calculates the Q value for two taxa.

   Receive: Distance matrix and two taxa.
   Return: Q value between the two taxa.
   '''

   return (len(distMat)-2)*calcDist(distMat, i, j) - \
           calcDistSum(distMat, i) - \
           calcDistSum(distMat, j)

def calcQMat(distMat):
   '''
   Calculates the entire Q matrix.

   Receive: Distance matrix.
   Return: Q matrix.
   '''

   q = numpy.zeros((len(distMat),len(distMat)), int)

   for i in range(1, len(distMat)):
      for j in range(i):
         q[i][j] = calcQ(distMat, i, j)

   return q

def calcDistOldNew(distMat, i, j):
   '''
   Calculates the distance from each of the old taxas (i,j) to the new combined taxa (ij).

   Receive: Distance matrix and the two old taxas.
   '''

   return (.5)*(calcDist(distMat, i, j)) + ((1./(2*(len(distMat)-2))) * \
          (calcDistSum(distMat,i) - calcDistSum(distMat, j)))

def calcDistNewOther(distMat, k, f, g):
   '''
   Calculates the distance between the new combined taxa (fg) and another taxa (k).

   Receive: The distance matrix, old taxas (f,g) and the other taxa (k).
   Return: Distance between (fg) and k.
   '''

   return (.5)*(calcDist(distMat,f,k) - calcDistOldNew(distMat, f, g)) + \
          (.5)*(calcDist(distMat,g,k) - calcDistOldNew(distMat, g, f))
   
def minQVal(q):
   '''
   Finds the minimum Q value, which will be the next combined taxa.

   Receive: Q matrix.
   Return: The minimum Q value and each of the taxa.
   '''

   iMin = 0
   jMin = 0
   qMin = 0

   for i in range(len(q)):
      for j in range(len(q)):
         if min(qMin, q[i][j]) == q[i][j]:
            qMin = q[i][j]
            iMin = i
            jMin = j

   if i > j:
      i, j = j, i

   return qMin, iMin, jMin

def doNeigJoin(mat, taxaList):
   '''
   Recursively executes the neighbor-join algorithm.  Ends when the size of the distance matrix is 1.

   Receive: The distance matrix, and a list of the taxa names.
   Return: If finished, the completed tree.  Otherwise, the next distance matrix and taxa list.

   Algorithm (when the matrix is greater than size 1):
      1. Calculate Q matrix.
      2. Find lowest value in Q matrix and combine corresponding taxa.
      3. Calculate new distance matrix that is one size smaller than previous step.
      4. Recursively call function until size is 1.
   '''

   if len(mat) == 1:
      return str(taxaList).replace(' ', '').replace('[','(').replace(']',')').replace('\'','')

   q = calcQMat(mat)

   minQ, taxaA, taxaB = minQVal(q)

   # initialize our new distance matrix
   newMat = numpy.zeros((len(mat)-1, len(mat)-1), float)

   # combine old taxa in taxa list to create new taxalist
   oldTaxaList = taxaList[:]
   oldTaxaList.remove(taxaList[taxaA])
   oldTaxaList.remove(taxaList[taxaB])
   newTaxaList = [[taxaList[taxaA], taxaList[taxaB]]] + oldTaxaList

   # calculate new distance matrix for new combined taxa values
   for i in range(1, len(newMat)):
      oldI = taxaList.index(newTaxaList[i])
      newMat[i][0] = calcDistNewOther(mat, oldI, taxaB, taxaA)

   # copy over everything else from old distance matrix
   for i in range(2, len(newMat)):
      for j in range(1, len(newMat)-1):
         oldI = taxaList.index(newTaxaList[i])
         oldJ = taxaList.index(newTaxaList[j])
         newMat[i][j] = mat[oldI][oldJ]

   return doNeigJoin(newMat, newTaxaList)

def getMaxInMatrix(mat):
   '''
   Finds the maximum value in a given matrix.

   Receive: A matrix.
   Recturn: The maximum value in the matrix.
   '''

   maxVal = 0

   for i in range(len(mat)):
      for j in range(len(mat)):
         if max(maxVal, mat[i][j]) == mat[i][j]:
            maxVal = mat[i][j]

   return maxVal

def normalizeMatrix(mat):
   '''
   Takes the BLOSUM alignment scores and normalizes the matrix so that the most similar sequences will
   have a smaller distance.  The most closely related taxa will have a score of 0.

   Receive: A matrix of BLOSUM alignment scores.
   Return: A distance matrix.
   '''

   maxVal = getMaxInMatrix(mat)

   for i in range(1, len(mat)):
      for j in range(i):
         mat[i][j] = maxVal - mat[i][j]

   return mat

def createTree(mat, taxaList):
   '''
   Calculates an unrooted tree based on a BLOSUM sequence alignment matrix.

   Receive: A matrix of BLOSUM scores and a taxa list.
   Return: A string representing the unrooted tree.
   '''

   mat = normalizeMatrix(mat)

   return doNeigJoin(mat, taxaList)

def main():
   '''
   The main function below is used for testing.
   '''

   taxaList = ['A', 'B', 'C', 'D', 'E', 'F']

   print taxaList

   mat = numpy.zeros((6,6), float)
   mat[1][0] = 5
   mat[2][0] = 4
   mat[2][1] = 7
   mat[3][0] = 7
   mat[3][1] = 10
   mat[3][2] = 7
   mat[4][0] = 6
   mat[4][1] = 9
   mat[4][2] = 6
   mat[4][3] = 5
   mat[5][0] = 8
   mat[5][1] = 11
   mat[5][2] = 8
   mat[5][3] = 9
   mat[5][4] = 8
  
   print doNeigJoin(mat, taxaList)

#main()