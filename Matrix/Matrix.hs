-----------------------------------------------------------------------------
-- |
-- Module      :  DSP.Matrix.Matrix
-- Copyright   :  (c) Matthew Donadio 2003
-- License     :  GPL
--
-- Maintainer  :  m.p.donadio@ieee.org
-- Stability   :  experimental
-- Portability :  portable
--
-- Basic matrix routines
--
-----------------------------------------------------------------------------

module Matrix.Matrix where

import Data.Array
import Data.Complex

-- | Matrix-matrix multiplication: A x B = C

mm_mult :: (Ix a, Integral a, Num b) => Array (a,a) b -- ^ A
	-> Array (a,a) b -- ^ B
	-> Array (a,a) b -- ^ C

mm_mult a b = if ac /= br
	      then error "mm_mult: inside dimensions inconsistent"
	      else array bnds [ ((i,j), mult i j) | (i,j) <- range bnds ]
    where mult i j = sum [ a!(i,k) * b!(k,j) | k <- [1..ac] ]
	  ((_,_),(ar,ac)) = bounds a
	  ((_,_),(br,bc)) = bounds b
	  bnds = ((1,1),(ar,bc))

-- | Matrix-vector multiplication: A x b = c

mv_mult :: (Ix a, Integral a, Num b) => Array (a,a) b -- ^ A
	-> Array a b -- ^ b
	-> Array a b -- ^ c

mv_mult a b = if ac /= br
	      then error "mv_mult: dimensions inconsistent"
	      else array bnds [ (i, mult i) | i <- range bnds ]
    where mult i = sum [ a!(i,k) * b!(k) | k <- [1..ac] ]
	  ((_,_),(ar,ac)) = bounds a
	  (_,br) = bounds b
	  bnds = (1,ar)

-- | Transpose of a matrix

m_trans :: (Ix a, Integral a, Num b) => Array (a,a) b -- ^ A
	-> Array (a,a) b -- ^ A^T

m_trans a = array bnds [ ((i,j), a!(j,i)) | (i,j) <- range bnds ]
    where (_,(m,n)) = bounds a
	  bnds = ((1,1),(n,m))

-- | Hermitian transpose (conjugate transpose) of a matrix

m_hermit :: (Ix a, Integral a, RealFloat b) => Array (a,a) (Complex b) -- ^ A
	 -> Array (a,a) (Complex b) -- ^ A^H

m_hermit a = array bnds [ ((i,j), conjugate (a!(j,i))) | (i,j) <- range bnds ]
    where (_,(m,n)) = bounds a
	  bnds = ((1,1),(n,m))
