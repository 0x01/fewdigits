{-# OPTIONS -XTypeOperators #-}
{-# OPTIONS -XTypeSynonymInstances #-}

module Data.Real.Integration
where

import Data.Real.Dyadic
import Data.Real.DReal
import Data.Real.Complete
import Data.Ratio
import Data.Multinomial

-- | half of (ceiling of 1 / 4-th root)
g :: Rational -> Int
g x | x == 0 = 0
    | True = (ceiling . (1/) . (2*) . toRational . sqrt . sqrt . fromRational) x
-- 'g' should be replaced by an exact function; the idea is to
-- approximate g to one digit and then check if this is indeed
-- the right value. If not it must be either g+1 or g-1
-- Alternatively compute as real real number and just add 2*eps
-- and take ceiling of that. This gives either g or g+1.

{-| In order to approximate an integral up to precision eps
  using Simpsons rule one must modify the interval at which
  the function is evaluated. This function computes the
  interval to the forth power: h^4 = 180 * epsilon / (M * (b-a))
  This is done using rationals and gives an exact answer,
  assuming b > a, m > 0, eps > 0.
-}
h4 :: Rational -> Rational -> Rational -> Gauge -> Rational
h4 a b m eps | (b-a) > 0 = 180 * eps / (m * (b-a) ^^ 5)
             | otherwise = 0

fromInt :: Num a => Int -> a
fromInt = fromInteger . toEnum

-- | To apply an uniform continuous function to a rational and
-- get a result with precision eps, we need to approximate the
-- input at precision mu eps
applyToRational :: (Dyadic :=> DReal) -> Rational -> DReal
applyToRational f = bind f . fromRational

simpson_steps :: Rational -> Rational -> Rational -> Gauge -> Int
simpson_steps a b m eps = 1 + 2 * (g $ h4 a b m eps)

dsimpson :: Rational -> Rational -> Rational -> (Dyadic :=> DReal) -> DReal
dsimpson a b m f eps = h3 * (sum $ map f' ps)
  where
    h3 :: Dyadic
    h3 = fromRational (h/3) eps
    ps = [(1,a)] ++  pts ++ [(1,b)]
    
    n :: Int
    n = 1 + 2 * (g $ h4 a b m eps)
    -- ^ number of points (a b excluded)
    h :: Rational
    h = (b - a) / (toInteger (n + 1) % 1)
    -- ^ distance between two points
    
    pts :: [(Int,Rational)]
    pts = take n $ zip (cycle [4,2]) ahh
    -- ^ actual points 4*x1, 2*x2, 4*x3, ...
    
    ahh = scanl (+) (a+h) (repeat h)

    f' (c, x) = (fromInt c) * (applyToRational f x e)
    e = eps / (toInteger n % 1)
    -- ^ we are applying f n times, so need precision eps/n        

-- example: dsimpson 1 2 1 (dLnUniformCts $ fromInteger 1 ) (1/10000000000000)
-- example: dsimpson 0 (31415729/5000000) 1 (dSinCts) (1/100)
-- example: dsimpson 0 1 1 idCts (1/1000) ≈ 1/2

subdivisions :: (Fractional a) => a -> a -> Int -> [a]
subdivisions lb ub n = (take n pts) ++ [ub]
  where pts = scanl (+) lb (repeat h)
        h = (ub - lb) / (fromInteger $ toEnum n)

polyapprox a b m f eps = lagrange (zip pts $ map f' pts)
  where pts = subdivisions a b (simpson_steps a b m eps)
        f' = toQ . ($ eps) . (bind f) . fromRational

symbolicSimpson a b m f eps = integrate $ polyapprox a b m f eps