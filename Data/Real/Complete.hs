{-# OPTIONS -XTypeOperators #-}

module Data.Real.Complete 
                (Gauge, Complete, unit, join, mapC, bind, bindCts, mapC2,
                 (:=>), mkUniformCts, modulus, forgetUniformCts,
                 (:=>>), mkContraction, lipschitzConst, forgetContraction,
                 asUniformCts, asContraction, constCts, idCts, o) where

import Control.Exception

infixr 9 :=>
infixr 9 :=>>

type Gauge = Rational 

-- Intended that d (f x) (f y) <= x + y
type Complete a = Gauge -> a

unit :: a -> Complete a
unit x eps = x

join :: (Complete (Complete a)) -> Complete a
join f eps = (f (eps/2)) (eps/2)

-- A uniformly continuous function on some subset of a to b
-- Hopefully the name of the function gives an indication of
-- the domain.
data a :=> b = UniformCts 
               {modulus :: (Gauge -> Gauge)
               ,forgetUniformCts :: (a -> b)}

mkUniformCts = UniformCts

mapC :: (a :=> b) -> Complete a -> Complete b
mapC (UniformCts mu f) x eps = f (x (mu eps))

bind :: (a :=> Complete b) -> Complete a -> Complete b
bind f x = join $ mapC f x

bindCts :: (a :=> Complete b) -> (Complete a :=> Complete b)
bindCts f = mkUniformCts (modulus f) (bind f)

mapC2 :: (a :=> b :=> c) -> Complete a -> Complete b -> Complete c
mapC2 f x y eps = (mapC approxf y) (eps/2)
 where
  approxf = (mapC f x) (eps/2)

o :: (b :=> c) -> (a :=> b) -> (a :=> c)
f `o` g = mkUniformCts mu h
 where
  mu = (modulus g) . (modulus f)
  h = (forgetUniformCts f) . (forgetUniformCts g)

constCts :: a -> b :=> a
constCts a = mkUniformCts (const undefined) (const a)

idCts :: a :=> Complete a
idCts = mkUniformCts id unit

data a :=>> b = ContractionMap
                { lipschitzConst :: Gauge
                , forgetContraction :: (a -> b) }

mkContraction c f = assert (c >= 0 && c < 1) $ ContractionMap c f

asUniformCts :: (a :=>> b) -> (a :=> b)
asUniformCts f = mkUniformCts (\e -> e / lipschitzConst f) (forgetContraction f)

asContraction :: (a :=> b) -> Gauge -> (a :=>> b)
asContraction f c = mkContraction c (forgetUniformCts f)