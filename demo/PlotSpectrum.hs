module Main where

import qualified Graphics.Gnuplot.Simple as Gnuplot

import qualified Numeric.Transform.Fourier.FFT as FT
import qualified Data.Array   as Array
import qualified Data.Complex as Complex

import qualified Sound.Sox.Read          as SoxRead
import qualified Sound.Sox.Format        as SoxFormat
import qualified Sound.Sox.Option.Format as SoxOption
import qualified Sound.Sox.Signal.List   as SoxList

import qualified Control.Monad.Exception.Asynchronous as AsyncExc


soundsDir :: FilePath
soundsDir =
   "samples/"

plotSpectrum :: [Double] -> IO ()
plotSpectrum x =
   let len = length x
   in  Gnuplot.plotList [] .
--       take 100 .
       take (div len 2) .
       Array.elems .
       fmap Complex.magnitude .
       FT.rfft .
       Array.listArray (0, len - 1) $
       x

main :: IO ()
main =
   plotSpectrum =<<
      fmap AsyncExc.result .
      SoxRead.withHandle2 SoxList.getContents =<<
      SoxRead.open (SoxOption.format SoxFormat.signedByte) (soundsDir++"BassDrum2")

{-
Result: BassDrum2 has two major peaks at periods 250 and 170 samples.
-}
