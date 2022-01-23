#!/usr/bin/python3

from typing import Match
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		# score = random.random()*100;
		# score = None
		# if banded:
		# 	score = self.edit_banded(seq1, seq2, align_length)[0]
		# else:
		score = self.edit(seq1, seq2, align_length)[0]
		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################					
		
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	def edit(self, seq1, seq2, align_length):
		#IF THE ALIGN_LENGTH IS LESS THAN THE SEQ THAN THAT IS THE LENGTH WE USE
		if align_length > len(seq1):
			len_seq1 = len(seq1)
		else:
			len_seq1 = align_length
		if align_length > len(seq2):
			len_seq2 = len(seq2)
		else:
			len_seq2 = align_length
		#define a 2D array that is |seq1| x |seq2|
		E = [[(None, None) for i in range(len_seq2+1)] for j in range(len_seq1+1)]
		#to change seq1 to the empty seq will take |seq1| deletes
		for i in range(0,len_seq1+1):
			if i == 0:
				E[i][0] = (i*INDEL, None)
			else:
				E[i][0] = (i*INDEL, 'top')
		#to change the empty seq into seq2 will take |seq2| inserts
		for j in range(0,len_seq2+1):
			if j == 0:
				E[0][j] = (j*INDEL, None)
			else:
				E[0][j] = (j*INDEL, 'side')

		for k in range(0,len_seq1):
			for m in range(0, len_seq2):
				diagonal = None
				side = None
				top = None
				if seq2[m] == seq1[k]:
					diagonal = MATCH + (E[k][m])[0]
				else:
					diagonal = SUB + (E[k][m])[0]
				side = (E[k+1][m])[0] + INDEL
				top = (E[k][m+1])[0] + INDEL
				E[k+1][m+1] = self.min3(diagonal, side, top)
		return E[len_seq1][len_seq2]
	
	def edit_banded(self, seq1, seq2, align_length):
		#IF THE ALIGN_LENGTH IS LESS THAN THE SEQ THAN THAT IS THE LENGTH WE USE
		if align_length > len(seq1):
			len_seq1 = len(seq1)
		else:
			len_seq1 = align_length
		if align_length > len(seq2):
			len_seq2 = len(seq2)
		else:
			len_seq2 = align_length
		#IF THE SEQUENCES HAVE A LEN DIFFERENCE > MAXINDELS THAN THERE IS NO POSSIBLE ALIGNMENT
		if abs(len_seq1-len_seq2) > MAXINDELS:
			return (math.inf, None)
		#CREATE AN nxlen_seq2 ARRAY CALLED E
		E = [[(None, None) for i in range(7)] for j in range(len_seq2+1)]

		for i in range(len_seq2+1):
			Val_for_empty_str = 0
			for j in range(7):
                #TOP EDGE CASES OF ARRAY E	
				if i == 0:
					if j < 3:
						E[i][j] = (math.inf, None)
					else:
						E[i][j] = (0 + Val_for_empty_str, 'Side')
						Val_for_empty_str+=INDEL
				elif i == 1 and j <= 2:
					if j < 2:
						E[i][j] = (math.inf, None)
					if j == 2:
						E[i][j] = (5, 'top')
				elif i == 2 and j <= 1:
					if j < 1:
						E[i][j] = (math.inf, None)
					if j == 1:
						E[i][j] = (10, 'top')
				elif i == 3 and j == 0:
						E[i][j] = (15, 'top')
				#BOTTOM EDGE CASES OF ARRAY E
				elif i == len_seq2 - 2 and j == 6:
						E[i][j] = (math.inf, None)
				elif i == len_seq2 -1 and j >= 5:
						E[i][j] = (math.inf, None)
				elif i == len_seq2 and j >= 4:
						E[i][j] = (math.inf, None)
				#REGULAR CASES
				else:
					diagonal = None
					side = None
					top = None
					char_seq = i - MAXINDELS + j - 1
					if seq2[char_seq] == seq1[char_seq]:
						diagonal = MATCH + (E[i-1][j])[0]
					else:
						diagonal = SUB + (E[i-1][j])[0]
					if j != 0:
						side = (E[i][j-1])[0] + INDEL
					if j + 1 <= 6:
						top = (E[i-1][j+1])[0] + INDEL
					if top == None:
						E[i][j] = self.min2_ds(diagonal, side)
					elif side == None:
						E[i][j] = self.min2_td(top, diagonal)
					else:
						E[i][j] = self.min3(diagonal, side, top)
		return E[len_seq2][0]

	def min3(self,diag,side,top):
		min = None
		direction = None
		if diag < side:
			min = diag
			direction = 'diagonal'
		else:
			min = side
			direction = 'side'
		if ((top <= min) and (direction =='diagonal')) or (top < min):
			min = top
			direction = 'top'
		return (min, direction)

	def min2_ds(self,diag,side):
		min = None
		direction = None
		if diag < side:
			min = diag
			direction = 'diagonal'
		else:
			min = side
			direction = 'side'
		return (min, direction)
	def min2_td(self, top, diagonal):
		min = None
		direction = None
		if diagonal < top:
			min = diagonal
			direction = 'diagonal'
		else:
			min = top
			direction = 'top'
		return (min, direction)

	def get_path(self, E):
		rows = len(E)
		cols = len(E[0])
		path = []
		path = self.get_align_helper(rows,cols,E)
		return path

	def get_path_helper(self, row, col, E):
		path=[]
		if row == 0 and col == 0:
			return path
		if (E[row][col])[1] == 'diagonal':
			path.append('diagonal')
			path.append(self.get_align_helper(row-1,col-1,E))
		elif (E[row][col])[1] == 'top':
			path.append('top')
			path.append(self.get_align_helper(row-1,col,E))
		elif (E[row][col])[1] == 'side':
			path.append('side')
			path.append(self.get_align_helper(row,col-1,E))
		return path

	def get_alignment(self, path, seq1, seq2):
		num_chars = len(path)	
		align_1 = ''
		align_2 = ''
		for i in range(num_chars):
			align_1 += '*'
			align_2 += '*'
		for i in range(num_chars-1, -1, -1):
			direction = path.pop(0)
			if direction == 'diagonal':
				align_1[i] = seq1[i]
				align_2[i] = seq2[i]
			if direction == 'top':
				align_1[i] = '-'
				align_2[i] = seq2[i]
			if direction == 'side':
				align_1 = seq1[i]
				align_2[i] = '-'
		return (align_1, align_2)
				


			# def edit_banded(self, seq1, seq2, align_length):
	# 	d = 3
	# 	band = 2*d+1
	# 	if align_length > len(seq1):
	# 		len_seq1 = len(seq1)
	# 	else:
	# 		len_seq1 = align_length
	# 	if align_length > len(seq2):
	# 		len_seq2 = len(seq2)
	# 	else:
	# 		len_seq2 = align_length
	# 	#define a 2D array that is |seq1| x |seq2|
	# 	# E =[[(None,None)]*len_seq1]*len_seq2
	# 	E = [[(None, None) for i in range(band)] for j in range(len_seq2+1)]
	# 	#to change seq1 to the empty seq will take |seq1| deletes
	# 	for i in range(0,5):
	# 		if i == 0:
	# 			E[i][0] = (i*INDEL, None)
	# 		else:
	# 			E[i][0] = (i*INDEL, 'top')
	# 	#to change the empty seq into seq2 will take |seq2| inserts
	# 	for j in range(0,5):
	# 		if j == 0:
	# 			E[0][j] = (j*INDEL, None)
	# 		else:
	# 			E[0][j] = (j*INDEL, 'side')

	# 	for i in range(1,5):
	# 		for j in range(1,4+i):
	# 			last = 3+i
	# 			if seq2[i] == seq1[j]:
	# 				diagonal = MATCH + (E[i][j])[0]
	# 			else:
	# 				diagonal = SUB + (E[i][j])[0]
	# 			side = (E[i][j+1])[0] + INDEL
	# 			if j != last:
	# 				top = (E[i+1][j])[0] + INDEL
	# 				E[i+1][j+1] = self.min(diagonal, side, top)
	# 			else:
	# 				E[i+1][j+1] = self.min(diagonal, side)

	# 	for k in range(0,len_seq2):
	# 		for m in range(0, band):
	# 			val_m = m
	# 			if k < 4:
	# 				m = 4+k
	# 			diagonal = None
	# 			side = None
	# 			top = None
	# 			if seq2[k] == seq1[m]:
	# 				diagonal = MATCH + (E[k][m])[0]
	# 			else:
	# 				diagonal = SUB + (E[k][m])[0]
	# 			side = (E[k+1][m])[0] + INDEL
	# 			top = (E[k][m+1])[0] + INDEL
	# 			E[k+1][m+1] = self.min(diagonal, side, top)
	# 	return E[len_seq2][len_seq1]	