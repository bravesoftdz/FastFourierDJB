unit ns_FFT_utils;

(* ********************************************************************
  *  ns_FFT_utils.pas                                                 *
  *                                                                   *
  *                                                                   *
  *  author    : Alexey V. Nikitayev                                  *
  *  email     : nikitayev@mail.ru                                    *
  *  web       : https://github.com/nikitayev                         *
  *  date      : 08-10-2015                                           *
  *                                                                   *
  *  The contents of this file are used with permission, subject to   *
  *  the Mozilla Public License Version 1.1 (the "License"); you may  *
  *  not use this file except in compliance with the License. You may *
  *  obtain a copy of the License at                                  *
  *  http://www.mozilla.org/MPL/MPL-1.1.html                          *
  *                                                                   *
  *  Software distributed under the License is distributed on an      *
  *  "AS IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or   *
  *  implied. See the License for the specific language governing     *
  *  rights and limitations under the License.                        *
  *                                                                   *
  *  Based on AlgLib and dspFastFourier                               *
  *  AlgLib available at http://www.alglib.net/download.php           *
  *  or: https://github.com/nikitayev/NeuroTests/                     *
  *                                                                   *
  *  dspFastFourier available at http://dsp-worx.de                   *
  *  or: https://github.com/nikitayev/FastFourierDJB/                 *
  *                                                                   *
  *  (C) 2015 Alexey V. Nikitaev <nikitayev@mail.ru>                  *
  *                                                                   *
  ******************************************************************** *)
  
interface

uses Math, dspConst, dspFastFourier, dspUtils, Ap, fft;

type
  { Used by TDCFFT to setup the FFT Size. }
  TDCFFTSize = (fts2, fts4, fts8, fts16, fts32, fts64, fts128, fts256, fts512, fts1024, fts2048, fts4096, fts8192,
    fts16384, fts32768);

function GetFFTLength(aFFTSize: TDCFFTSize): Integer;

procedure EKG2FFT_LF(const aData: TReal1DArray; aFFTSize: TDCFFTSize; var aFFT: TComplex1DArray;
  aWindowMode: TWindowMode);
procedure EKG2FFT_LF_AlgLib(const aData: TReal1DArray; aFFTSize: TDCFFTSize; var aFFT: TComplex1DArray;
  aWindowMode: TWindowMode);
procedure FFT2EKG(const aFFT: TComplex1DArray; aFFTSize: TDCFFTSize; var aData: TComplex1DArray);
procedure FFT2EKG_AlgLib(const aFFT: TComplex1DArray; aFFTSize: TDCFFTSize; var aData: TComplex1DArray);

procedure FFT_CalcFreqAndPhase(aFFTSize: TDCFFTSize; var arr: TComplex1DArray; nSamplesPerSec: AlglibFloat;
  out freq, amp, phase: TReal1DArray);

procedure FFT_CalcAmpAndFreq(aFFTSize: TDCFFTSize; var arr: TComplex1DArray; aTotalTimeSec: AlglibFloat;
  out freq, amp: TReal1DArray);

function FFTDeleteLow(var aFFT: TComplex1DArray; aLowLimit: AlglibFloat): Integer;
procedure FFTPack(var aFFT: TComplex1DArray; aFFTSize: TDCFFTSize);
procedure FFTUnPack(var aFFT: TComplex1DArray; aFFTSize: TDCFFTSize);

procedure EKG2FFT_HF(const aEKG: TReal1DArray; aFFTSize: TDCFFTSize; var aFFT: TReal1DArray);

implementation

function GetFFTLength(aFFTSize: TDCFFTSize): Integer;
begin
  result := 1 shl (Integer(aFFTSize) + 1);
end;

function FFTBuff_AllocMem(aFFTSize: TDCFFTSize): PComplexArray;
begin
  result := AllocMem(GetFFTLength(aFFTSize) * SizeOf(TComplex));
end;

procedure FFTApplyWindow(var aFFT: TComplex1DArray; aWindowMode: TWindowMode);
var
  i: Integer;
begin
  for i := 0 to High(aFFT) do
    aFFT[i].X := GetWindowingValue(aFFT[i].X, i, Length(aFFT), aWindowMode);
end;

// вычисление FFT с помощью библиотеки dspFastFourier
procedure EKG2FFT_LF(const aData: TReal1DArray; aFFTSize: TDCFFTSize; var aFFT: TComplex1DArray;
  aWindowMode: TWindowMode);
var
  i: Integer;
  zFFTLength: Integer;
  zFFT: PComplexArray;
  zMaxValue, zThresholdValue: AlglibFloat;
begin
  zFFTLength := GetFFTLength(aFFTSize);
  zFFT := FFTBuff_AllocMem(aFFTSize);
  SetLength(aFFT, zFFTLength);

  try
    // занулим
    for i := 0 to zFFTLength - 1 do
    begin
      zFFT[i].re := aData[i];
      zFFT[i].im := 0;
    end;

    // применим окно выбранного типа
    for i := 0 to zFFTLength - 1 do
      zFFT[i].re := GetWindowingValue(zFFT[i].re, i, zFFTLength, aWindowMode);

    dspDoFFT(zFFT, zFFTLength, False, True, True);

    for i := 0 to zFFTLength - 1 do
    begin
      aFFT[i].X := zFFT[i].re;
      aFFT[i].Y := zFFT[i].im;
    end;

  finally
    FreeMem(zFFT);
  end;
end;

// усреднение исходных данных - дл€ возможности оценки сверхнизких частот
procedure EKG2FFT_LF_AlgLib(const aData: TReal1DArray; aFFTSize: TDCFFTSize; var aFFT: TComplex1DArray;
  aWindowMode: TWindowMode);
var
  i: Integer;
  zFFTLength: Integer;
begin
  zFFTLength := GetFFTLength(aFFTSize);
  SetLength(aFFT, zFFTLength);

  for i := 0 to zFFTLength - 1 do
  begin
    aFFT[i].X := aData[i];
    aFFT[i].Y := 0;
  end;

  FFTApplyWindow(aFFT, aWindowMode);

  FFTC1D(aFFT, zFFTLength);
end;

procedure FFT2EKG(const aFFT: TComplex1DArray; aFFTSize: TDCFFTSize; var aData: TComplex1DArray);
var
  i: Integer;
  zFFT: PComplexArray;
  zFFTLength: Integer;
  zMaxValue, zThresholdValue: AlglibFloat;
begin
  zFFTLength := GetFFTLength(aFFTSize);
  zFFT := FFTBuff_AllocMem(aFFTSize);
  SetLength(aData, zFFTLength);

  try

    for i := 0 to High(aFFT) do
    begin
      zFFT[i].re := aFFT[i].X;
      zFFT[i].im := aFFT[i].Y;
    end;

    dspDoFFT(zFFT, zFFTLength, True, True, True);

    for i := 0 to High(aFFT) do
    begin
      aData[i].X := zFFT[i].re;
      aData[i].Y := zFFT[i].im;
    end;
  finally
    FreeMem(zFFT);
  end;
end;

procedure FFT2EKG_AlgLib(const aFFT: TComplex1DArray; aFFTSize: TDCFFTSize; var aData: TComplex1DArray);
var
  i: Integer;
  zFFTLength: Integer;
  zMaxValue, zThresholdValue: AlglibFloat;
begin
  zFFTLength := GetFFTLength(aFFTSize);
  SetLength(aData, zFFTLength);

  for i := 0 to High(aFFT) do
  begin
    aData[i] := aFFT[i];
  end;

  FFTC1DInv(aData, zFFTLength);
end;

// вз€то тут: http://psi-logic.narod.ru/fft/fftg.htm
// arr -  Ёто комплексные числа, результат пр€мого дискретного преобразовани€ ‘урье
procedure FFT_CalcFreqAndPhase(aFFTSize: TDCFFTSize; var arr: TComplex1DArray; nSamplesPerSec: AlglibFloat;
  out freq, amp, phase: TReal1DArray);
var
  N: Integer;

  // это индекс комплексного числа в массиве arr
  i: Integer;
  Nmax: Integer;

  // это частота дискретизации
  abs2: AlglibFloat;
  re, im: AlglibFloat;

begin
  N := GetFFTLength(aFFTSize);

  SetLength(arr, N);

  // убираем зеркальный эффект, просто отбрасыва€ вторую половину
  Nmax := (N + 1) div 2;

  // мы хотим получить массив гармоник.
  // это массив амплитуд, массив частот и массив фаз дл€ каждой гармоники
  SetLength(freq, Nmax);
  SetLength(amp, Nmax);
  SetLength(phase, Nmax);

  // получаем остальные гармоники
  for i := 0 to Nmax - 1 do
  begin
    re := arr[i].X;
    im := arr[i].Y;

    // это квадрат модул€ комплексного числа arr[i]
    abs2 := re * re + im * im;

    // вычисл€ем апмлитуду. 2.0 - дл€ устранени€ зеркального эффекта
    amp[i] := 2.0 * sqrt(abs2) / N;

    // вычисл€ем фазу косинуса в радианах
    phase[i] := ArcTan2(im, re);

    // преобразуем косинус в синус. M_PI2 = пи/2, M_PI = пи
    // в результате фаза будет в диапазоне от -пи/2 до +пи/2
    phase[i] := phase[i] + (PI * 0.5);
    if (phase[i] > PI) then
      phase[i] := phase[i] - 2 * PI;

    // можно еще преобразовать радианы в градусы
    phase[i] := phase[i] * (180.0 / PI);

    // получаем частоту
    freq[i] := (nSamplesPerSec * i) / N;
  end;
end;

procedure FFT_CalcAmpAndFreq(aFFTSize: TDCFFTSize; var arr: TComplex1DArray; aTotalTimeSec: AlglibFloat;
  out freq, amp: TReal1DArray);
var
  N: Integer;

  // это индекс комплексного числа в массиве arr
  i: Integer;
  Nmax: Integer;

  // это частота дискретизации
  abs2: AlglibFloat;
  re, im: AlglibFloat;

begin
  N := GetFFTLength(aFFTSize);

  SetLength(arr, N);

  // убираем зеркальный эффект, просто отбрасыва€ вторую половину
  Nmax := (N + 1) div 2;

  // мы хотим получить массив гармоник.
  // это массив амплитуд, массив частот и массив фаз дл€ каждой гармоники
  SetLength(freq, Nmax);
  SetLength(amp, Nmax);

  // получаем остальные гармоники
  for i := 0 to Nmax - 1 do
  begin
    re := arr[i].X;
    im := arr[i].Y;

    // это квадрат модул€ комплексного числа arr[i]
    abs2 := re * re + im * im;

    // вычисл€ем апмлитуду. 2.0 - дл€ устранени€ зеркального эффекта
    amp[i] := 2.0 * sqrt(abs2) / N;

    // получаем частоту
    freq[i] := i / aTotalTimeSec;
  end;
end;

// уберЄм значени€, которые меньше заданных
function FFTDeleteLow(var aFFT: TComplex1DArray; aLowLimit: AlglibFloat): Integer;
var
  i: Integer;
begin
  result := 0;
  aLowLimit := sqr(aLowLimit) * sqr(Length(aFFT));
  for i := 1 to High(aFFT) do
    with aFFT[i] do
      if (X * X + Y * Y < aLowLimit) then
      begin
        X := 0;
        Y := 0;
        inc(result);
      end;
end;

// упаковка результатов пр€мого преобразовани€
procedure FFTPack(var aFFT: TComplex1DArray; aFFTSize: TDCFFTSize);
begin
  SetLength(aFFT, (GetFFTLength(aFFTSize) div 2) + 1);
end;

// распаковка результатов пр€мого преобразовани€
procedure FFTUnPack(var aFFT: TComplex1DArray; aFFTSize: TDCFFTSize);
var
  i: Integer;
  N2: Integer;
  function InvertIm(const aValue: Complex): Complex; inline;
  begin
    result.X := aValue.X;
    result.Y := -aValue.Y;
  end;

begin
  N2 := GetFFTLength(aFFTSize) div 2;
  SetLength(aFFT, GetFFTLength(aFFTSize));
  for i := 1 to N2 - 1 do
  begin
    aFFT[N2 + i] := InvertIm(aFFT[N2 - i]);
  end;
end;

// получение усреднЄнного спектра частот дл€ интервалов
procedure EKG2FFT_HF(const aEKG: TReal1DArray; aFFTSize: TDCFFTSize; var aFFT: TReal1DArray);
var
  i, j, zCount, zOffset: Integer;
  zFFTLength: Integer;
  zFFT, zFFTSwap: PComplexArray;
  zK: AlglibFloat;
begin
  zFFTLength := GetFFTLength(aFFTSize);
  zFFT := FFTBuff_AllocMem(aFFTSize);
  zFFTSwap := FFTBuff_AllocMem(aFFTSize);
  SetLength(aFFT, zFFTLength div 2);

  try
    zCount := Length(aEKG) div zFFTLength;

    for j := 0 to zCount - 1 do
    begin
      zOffset := j * zFFTLength;
      for i := 0 to zFFTLength - 1 do
      begin
        zFFT[i].re := aEKG[i + zOffset];
        zFFT[i].im := 0;
      end;

      for i := 0 to zFFTLength - 1 do
        zFFT[i].re := GetWindowingValue(zFFT[i].re, i, zFFTLength, TWindowMode.wmHamming);

      dspDoFFT(zFFT, zFFTLength, False, True, True, zFFTSwap);

      for i := 0 to (zFFTLength shr 1) - 1 do
      begin
        aFFT[i] := aFFT[i] + FFTSum(zFFT[i].re, zFFT[i].im);
      end;
    end;

    // найдЄм среднее
    zK := 1.0 / zCount;
    for i := 0 to (zFFTLength shr 1) - 1 do
    begin
      aFFT[i] := aFFT[i] * zK;
    end;
    aFFT[0] := 0;
  finally
    FreeMem(zFFT);
    FreeMem(zFFTSwap);
  end;
end;

end.
