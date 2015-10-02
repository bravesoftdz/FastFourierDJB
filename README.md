# FastFourierDJB
Modification of FastFourierDJB
                                                               
 author    : Milenko Mitrovic                                    
 email     : dcoder@dsp-worx.de                                  
 web       : http://dsp-worx.de                                  
 date      : 24-07-2003                                          
                                                                 
 Based on D.J.Bernsteinґs Split Radix FFT v0.76.                 
 http://cr.yp.to/djbfft.html                                     
                                                                 
 The contents of this file are used with permission, subject to  
 the Mozilla Public License Version 1.1 (the "License"); you may 
 not use this file except in compliance with the License. You may
 obtain a copy of the License at                                 
 http://www.mozilla.org/MPL/MPL-1.1.html                         
                                                                 
 Software distributed under the License is distributed on an     
 "AS IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or  
 implied. See the License for the specific language governing    
 rights and limitations under the License.                       
                                                                 
 (C) 2003, 2004 Milenko Mitrovic <dcoder@dsp-worx.de>			 
																 
 Modifications by Alexey V. Nikitaev							 
 email     : nikitayev@mail.ru 									 
 web       : https://github.com/nikitayev/						 
 date      : 02-10-2015 										 
																 
 Sample:														 
																 
// применим окно Хэмминга										 
for i := 0 to zFFTLength - 1 do									 
  zFFT[i].re := GetWindowingValue(zFFT[i].re, i, zFFTLength, TWindowMode.wmHamming);
																 
dspDoFFT(zFFT, zFFTLength, False, True, True, zFFTSwap);		 
																 
zMaxValue := 0;													 
for i := 0 to (zFFTLength shr 1) - 1 do							 
begin															 
  aFFT[i] := FFTSum(zFFT[i].re, zFFT[i].im);					 
end;															 
