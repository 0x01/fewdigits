{-# OPTIONS -XTypeOperators #-}
{-# OPTIONS -XTypeSynonymInstances #-}

module Data.Real.Picard
where

import Data.Real.Dyadic
import Data.Real.DReal
import Data.Real.Integration
import Data.Real.Complete
import Data.Ratio

-- | compose a function n times
ncomp :: Int -> (a -> a) -> (a -> a)
ncomp 1 = id
ncomp n = \f -> f . ncomp (n - 1) f

-- contraction operator on the space of uniform cts real fns
type Operator a b = (a :=> b) :=>> (a :=> b)

-- | computes for a given contraction and precision the required
-- picard steps, ie, nr of compositions in (F o F o ... o F)
picard_steps :: Operator a b -> Gauge -> Int
picard_steps f eps = ceiling $ logBase c e -- TODO use precise log!
  where c = fromRational . lipschitzConst $ f 
        e = fromRational eps

picard :: Operator a (Complete a) -> a -> Complete a
picard op x eps = (forgetUniformCts opn) x eps
   where n = picard_steps op eps
         opn = ncomp n (forgetContraction op) idCts

-- | The initial value problem y'=y for y(0)=1 as intergral eqn
-- The exp function is the obvious solution to this IVP.
expOp :: Operator Dyadic DReal
expOp = mkContraction (1/2) fCts
  where fCts u = mkUniformCts (/2) (\t -> 1 + dsimpson 0 (toQ t) 3 u)

-- example: compute (exp .5): picard expOp (Dyadic 1 (- 1)) 1/10

{- bla
symbolicPicard :: Fractional a => Operator a (Complete a) -> a -> Polynomial a
symbolicPicard op x eps =  foldl (flip (.)) toP (replicate n integrate)
  where toP :: a -> Polynomial a
        toP x = symbolicSimpson 0 x 1 (asUniformCts op) eps
        n = max 0 $ picard_steps op eps - 1
-}