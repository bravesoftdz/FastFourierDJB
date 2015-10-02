(* ********************************************************************
  *  FastFourierDJB.pas                                               *
  *                                                                   *
  *  author    : Milenko Mitrovic                                     *
  *  email     : dcoder@dsp-worx.de                                   *
  *  web       : http://dsp-worx.de                                   *
  *  date      : 24-07-2003                                           *
  *                                                                   *
  *  Based on D.J.Bernstein´s Split Radix FFT v0.76.                  *
  *  http://cr.yp.to/djbfft.html                                      *
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
  *  (C) 2003, 2004 Milenko Mitrovic <dcoder@dsp-worx.de>             *
  *                                                                   *
  ******************************************************************** *)

{$I dspFFTConfig.inc}

unit FastFourierDJB;

interface

uses dspConst, Math;

var
  RevBits2: array [0 .. 1] of WORD;
  RevBits4: array [0 .. 3] of WORD;
  RevBits8: array [0 .. 7] of WORD;
  RevBits16: array [0 .. 15] of WORD;
  RevBits32: array [0 .. 31] of WORD;
  RevBits64: array [0 .. 63] of WORD;
  RevBits128: array [0 .. 127] of WORD;
  RevBits256: array [0 .. 255] of WORD;
  RevBits512: array [0 .. 511] of WORD;
  RevBits1024: array [0 .. 1023] of WORD;
  RevBits2048: array [0 .. 2047] of WORD;
  RevBits4096: array [0 .. 4095] of WORD;
  RevBits8192: array [0 .. 8191] of WORD;
  RevBits16384: array [0 .. 16383] of WORD;

var
  d16: array [0 .. 2] of TComplex;
  d32: array [0 .. 6] of TComplex;
  d64: array [0 .. 14] of TComplex;
  d128: array [0 .. 30] of TComplex;
  d256: array [0 .. 62] of TComplex;
  d512: array [0 .. 126] of TComplex;
  d1024: array [0 .. 126] of TComplex;
  d2048: array [0 .. 254] of TComplex;
  d4096: array [0 .. 510] of TComplex;
  d8192: array [0 .. 1022] of TComplex;
  d16384: array [0 .. 2046] of TComplex;

const
  SqrtHalf: Float = 0.70710678118654752440084436210484903;

  scalec2: Float = 0.5;
  scalec4: Float = 0.25;
  scalec8: Float = 0.125;
  scalec16: Float = 0.0625;
  scalec32: Float = 0.03125;
  scalec64: Float = 0.015625;
  scalec128: Float = 0.0078125;
  scalec256: Float = 0.00390625;
  scalec512: Float = 0.001953125;
  scalec1024: Float = 0.0009765625;
  scalec2048: Float = 0.00048828125;
  scalec4096: Float = 0.000244140625;
  scalec8192: Float = 0.0001220703125;
  scalec16384: Float = 0.00006103515625;

procedure fft2(a: PComplexArray);
procedure fft4(a: PComplexArray);
procedure fft8(a: PComplexArray);
procedure fft16(a: PComplexArray);
procedure fft32(a: PComplexArray);
procedure fft64(a: PComplexArray);
procedure fft128(a: PComplexArray);
procedure fft256(a: PComplexArray);
procedure fft512(a: PComplexArray);
procedure fft1024(a: PComplexArray);
procedure fft2048(a: PComplexArray);
procedure fft4096(a: PComplexArray);
procedure fft8192(a: PComplexArray);
procedure fft16384(a: PComplexArray);

procedure ifft2(a: PComplexArray);
procedure ifft4(a: PComplexArray);
procedure ifft8(a: PComplexArray);
procedure ifft16(a: PComplexArray);
procedure ifft32(a: PComplexArray);
procedure ifft64(a: PComplexArray);
procedure ifft128(a: PComplexArray);
procedure ifft256(a: PComplexArray);
procedure ifft512(a: PComplexArray);
procedure ifft1024(a: PComplexArray);
procedure ifft2048(a: PComplexArray);
procedure ifft4096(a: PComplexArray);
procedure ifft8192(a: PComplexArray);
procedure ifft16384(a: PComplexArray);

procedure scale2(a: PComplexArray);
procedure scale4(a: PComplexArray);
procedure scale8(a: PComplexArray);
procedure scale16(a: PComplexArray);
procedure scale32(a: PComplexArray);
procedure scale64(a: PComplexArray);
procedure scale128(a: PComplexArray);
procedure scale256(a: PComplexArray);
procedure scale512(a: PComplexArray);
procedure scale1024(a: PComplexArray);
procedure scale2048(a: PComplexArray);
procedure scale4096(a: PComplexArray);
procedure scale8192(a: PComplexArray);
procedure scale16384(a: PComplexArray);

procedure mul2(a, b: PComplexArray);
procedure mul4(a, b: PComplexArray);
procedure mul8(a, b: PComplexArray);
procedure mul16(a, b: PComplexArray);
procedure mul32(a, b: PComplexArray);
procedure mul64(a, b: PComplexArray);
procedure mul128(a, b: PComplexArray);
procedure mul256(a, b: PComplexArray);
procedure mul512(a, b: PComplexArray);
procedure mul1024(a, b: PComplexArray);
procedure mul2048(a, b: PComplexArray);
procedure mul4096(a, b: PComplexArray);
procedure mul8192(a, b: PComplexArray);
procedure mul16384(a, b: PComplexArray);

implementation

procedure TRANSFORM(a0, a1, a2, a3: PComplex; wre, wim: Float);
var
  t1, t2, t3, t4, t5, t6: Float;
begin
  t5 := a0.re - a2.re;
  t6 := a0.im - a2.im;
  t4 := a1.re - a3.re;
  t3 := a1.im - a3.im;

  a0.re := a0.re + a2.re;
  a0.im := a0.im + a2.im;
  a1.im := a1.im + a3.im;
  a1.re := a1.re + a3.re;

  t1 := t5 - t3;
  t2 := t6 - t4;
  t5 := t5 + t3;
  t6 := t6 + t4;

  a2.re := (t1 * wre) - (t6 * wim);
  a3.im := (t2 * wre) - (t5 * wim);
  a2.im := (t6 * wre) + (t1 * wim);
  a3.re := (t5 * wre) + (t2 * wim);
end;

procedure TRANSFORMHALF(a0, a1, a2, a3: PComplex);
var
  t1, t2, t3, t4, t5: Float;
begin
  t5 := a0.re - a2.re;
  t3 := a0.im - a2.im;
  t1 := a1.re - a3.re;
  t4 := a1.im - a3.im;

  a0.re := a0.re + a2.re;
  a0.im := a0.im + a2.im;
  a1.re := a1.re + a3.re;
  a1.im := a1.im + a3.im;

  t2 := t3 + t1;
  t3 := t3 - t1;
  t1 := t5 - t4;
  t5 := t5 + t4;

  a3.re := (t5 + t3) * SqrtHalf;
  a3.im := (t3 - t5) * SqrtHalf;
  a2.re := (t1 - t2) * SqrtHalf;
  a2.im := (t2 + t1) * SqrtHalf;
end;

procedure TRANSFORMZERO(a0, a1, a2, a3: PComplex);
var
  t1, t2, t3, t4: Float;
begin
  t1 := a0.re - a2.re;
  t2 := a0.im - a2.im;
  t3 := a1.re - a3.re;
  t4 := a1.im - a3.im;

  a0.re := a0.re + a2.re;
  a0.im := a0.im + a2.im;
  a1.re := a1.re + a3.re;
  a1.im := a1.im + a3.im;

  a2.re := t1 - t4;
  a3.im := t2 - t3;
  a3.re := t1 + t4;
  a2.im := t2 + t3;
end;

procedure c2(a: PComplexArray);
var
  t1, t2, t3, t4: Float;
begin
  t1 := a^[0].re + a^[1].re;
  t2 := a^[0].im + a^[1].im;
  t3 := a^[0].re - a^[1].re;
  t4 := a^[0].im - a^[1].im;

  a^[0].re := t1;
  a^[0].im := t2;
  a^[1].re := t3;
  a^[1].im := t4;
end;

procedure c4(a: PComplexArray);
var
  t1, t2, t3, t4, t5, t6, t7, t8: Float;
begin
  t3 := a^[0].re + a^[2].re;
  t2 := a^[0].im - a^[2].im;
  t4 := a^[1].re + a^[3].re;
  t6 := a^[1].im + a^[3].im;
  t1 := a^[0].re - a^[2].re;
  t7 := a^[0].im + a^[2].im;
  t5 := a^[1].re - a^[3].re;
  t8 := a^[1].im - a^[3].im;

  a^[0].re := t3 + t4;
  a^[0].im := t7 + t6;
  a^[1].re := t3 - t4;
  a^[1].im := t7 - t6;
  a^[2].im := t2 + t5;
  a^[2].re := t1 - t8;
  a^[3].im := t2 - t5;
  a^[3].re := t1 + t8;
end;

procedure c8(a: PComplexArray);
var
  t1, t2, t3, t4, t5, t6: Float;
begin
  t1 := a^[0].re - a^[4].re;
  t2 := a^[0].im - a^[4].im;
  t3 := a^[2].re - a^[6].re;
  t4 := a^[2].im - a^[6].im;

  a^[0].re := a^[0].re + a^[4].re;
  a^[0].im := a^[0].im + a^[4].im;
  a^[2].re := a^[2].re + a^[6].re;
  a^[2].im := a^[2].im + a^[6].im;

  a^[4].re := t1 - t4;
  a^[4].im := t2 + t3;
  a^[6].re := t1 + t4;
  a^[6].im := t2 - t3;

  t3 := a^[1].re - a^[5].re;
  t4 := a^[1].im - a^[5].im;
  t5 := a^[3].re - a^[7].re;
  t6 := a^[3].im - a^[7].im;

  a^[1].re := a^[1].re + a^[5].re;
  a^[1].im := a^[1].im + a^[5].im;
  a^[3].re := a^[3].re + a^[7].re;
  a^[3].im := a^[3].im + a^[7].im;

  t1 := t3 - t6;
  t3 := t3 + t6;
  t2 := t4 - t5;
  t4 := t4 + t5;
  t6 := (t1 - t4) * SqrtHalf;
  t1 := (t1 + t4) * SqrtHalf;
  t5 := (t2 - t3) * SqrtHalf;
  t2 := (t2 + t3) * SqrtHalf;

  a^[5].re := a^[4].re - t6;
  a^[5].im := a^[4].im - t1;
  a^[4].re := a^[4].re + t6;
  a^[4].im := a^[4].im + t1;

  a^[7].re := a^[6].re - t2;
  a^[7].im := a^[6].im - t5;
  a^[6].re := a^[6].re + t2;
  a^[6].im := a^[6].im + t5;
  c4(a);
end;

procedure c16(a: PComplexArray);
begin
  TRANSFORMZERO(@a^[0], @a^[4], @a^[8], @a^[12]);
  TRANSFORM(@a^[1], @a^[5], @a^[9], @a^[13], d16[0].re, d16[0].im);
  TRANSFORMHALF(@a^[2], @a^[6], @a^[10], @a^[14]);
  TRANSFORM(@a^[3], @a^[7], @a^[11], @a^[15], d16[0].im, d16[0].re);
  c4(@a^[8]);
  c4(@a^[12]);
  c8(a);
end;

procedure cpass(a, w: PComplexArray; n: integer);
var
  k: integer;
  z: integer;
begin
  TRANSFORMZERO(@a^[0], @a^[2 * n], @a^[4 * n], @a^[6 * n]);
  k := 2 * n - 1;
  z := 0;
  while k > 0 do
  begin
    TRANSFORM(@a^[1 + z], @a^[2 * n + 1 + z], @a^[4 * n + 1 + z], @a^[6 * n + 1 + z], w[z].re, w^[z].im);
    inc(z);
    dec(k);
  end;
end;

procedure c32(a: PComplexArray);
begin
  cpass(a, @d32, 4);
  c8(@a^[16]);
  c8(@a^[24]);
  c16(a);
end;

procedure c64(a: PComplexArray);
begin
  cpass(a, @d64, 8);
  c16(@a^[32]);
  c16(@a^[48]);
  c32(a);
end;

procedure c128(a: PComplexArray);
begin
  cpass(a, @d128, 16);
  c32(@a^[64]);
  c32(@a^[96]);
  c64(a);
end;

procedure c256(a: PComplexArray);
begin
  cpass(a, @d256, 32);
  c64(@a^[128]);
  c64(@a^[192]);
  c128(a);
end;

procedure c512(a: PComplexArray);
begin
  cpass(a, @d512, 64);
  c128(@a^[384]);
  c128(@a^[256]);
  c256(a);
end;

procedure cpassbig(a, w: PComplexArray; n: integer);
var
  k: integer;
  z: integer;
begin
  TRANSFORMZERO(@a^[0], @a^[2 * n], @a^[4 * n], @a^[6 * n]);

  k := n - 1;

  z := 0;
  while k > 0 do
  begin
    TRANSFORM(@a^[1 + z], @a^[2 * n + 1 + z], @a^[4 * n + 1 + z], @a^[6 * n + 1 + z], w[z].re, w^[z].im);
    inc(z);
    dec(k);
  end;

  TRANSFORMHALF(@a^[1 + z], @a^[2 * n + 1 + z], @a^[4 * n + 1 + z], @a^[6 * n + 1 + z]);
  inc(z);
  k := n - 1;

  while k > 0 do
  begin
    TRANSFORM(@a^[1 + z], @a^[2 * n + 1 + z], @a^[4 * n + 1 + z], @a^[6 * n + 1 + z], w[k - 1].im, w^[k - 1].re);
    inc(z);
    dec(k);
  end;
end;

procedure c1024(a: PComplexArray);
begin
  cpassbig(a, @d1024, 128);
  c256(@a^[768]);
  c256(@a^[512]);
  c512(a);
end;

procedure c2048(a: PComplexArray);
begin
  cpassbig(a, @d2048, 256);
  c512(@a^[1536]);
  c512(@a^[1024]);
  c1024(a);
end;

procedure c4096(a: PComplexArray);
begin
  cpassbig(a, @d4096, 512);
  c1024(@a^[3072]);
  c1024(@a^[2048]);
  c2048(a);
end;

procedure c8192(a: PComplexArray);
begin
  cpassbig(a, @d8192, 1024);
  c2048(@a^[6144]);
  c2048(@a^[4096]);
  c4096(a);
end;

procedure c16384(a: PComplexArray); inline;
begin
  cpassbig(a, @d16384, 2048);
  c4096(@a^[12288]);
  c4096(@a^[8192]);
  c4096(a);
end;

procedure UNTRANSFORM(a0, a1, a2, a3: PComplex; wre, wim: Float);
var
  t1, t2, t3, t4, t5, t6: Float;
begin
  t1 := a2.re * wre;
  t3 := a2.im * wim;
  t5 := a3.re * wre;
  t2 := a2.im * wre;
  t4 := a2.re * wim;
  t6 := a3.im * wre;
  t1 := t1 + t3;
  t5 := t5 - a3.im * wim;
  t2 := t2 - t4;
  t6 := t6 + a3.re * wim;
  t3 := t5 + t1;
  t5 := t5 - t1;
  t4 := t2 - t6;
  t6 := t6 + t2;
  a2.im := a0.im - t6;
  a3.re := a1.re - t4;
  a2.re := a0.re - t3;
  a1.re := t4 + a1.re;
  a0.im := t6 + a0.im;
  a0.re := t3 + a0.re;
  a3.im := a1.im - t5;
  a1.im := t5 + a1.im;
end;

procedure UNTRANSFORMHALF(a0, a1, a2, a3: PComplex);
var
  t1, t2, t3, t4, t5, t6: Float;
begin
  t1 := (a2.re + a2.im) * SqrtHalf;
  t2 := (a2.im - a2.re) * SqrtHalf;
  t3 := (a3.re - a3.im) * SqrtHalf;
  t4 := (a3.im + a3.re) * SqrtHalf;
  t6 := t3 - t1;
  t5 := t2 - t4;
  t1 := t1 + t3;
  t2 := t2 + t4;
  a3.im := a1.im - t6;
  a1.im := t6 + a1.im;
  a3.re := a1.re - t5;
  a1.re := t5 + a1.re;
  a2.re := a0.re - t1;
  a0.re := t1 + a0.re;
  a2.im := a0.im - t2;
  a0.im := t2 + a0.im;
end;

procedure UNTRANSFORMZERO(a0, a1, a2, a3: PComplex);
var
  t1, t2, t3, t4: Float;
begin
  t1 := a2.re + a3.re;
  t2 := a2.im + a3.im;
  t3 := a2.im - a3.im;
  t4 := a3.re - a2.re;
  a2.re := a0.re - t1;
  a2.im := a0.im - t2;
  a3.re := a1.re - t3;
  a3.im := a1.im - t4;
  a0.re := a0.re + t1;
  a0.im := a0.im + t2;
  a1.re := a1.re + t3;
  a1.im := a1.im + t4;
end;

procedure u4(a: PComplexArray);
var
  t1, t2, t3, t4, t5, t6, t7, t8: Float;
begin
  t1 := a^[0].re + a^[1].re;
  t2 := a^[0].im + a^[1].im;
  t3 := a^[0].re - a^[1].re;
  t4 := a^[0].im - a^[1].im;
  t6 := a^[3].re + a^[2].re;
  t8 := a^[3].im + a^[2].im;
  t5 := a^[3].re - a^[2].re;
  t7 := a^[2].im - a^[3].im;

  a^[0].re := t1 + t6;
  a^[2].re := t1 - t6;
  a^[0].im := t2 + t8;
  a^[1].re := t3 + t7;
  a^[3].re := t3 - t7;
  a^[1].im := t4 + t5;
  a^[3].im := t4 - t5;
  a^[2].im := t2 - t8;
end;

procedure u8(a: PComplexArray);
var
  t1, t2, t3, t4, t5, t6: Float;
begin
  u4(a);

  t1 := a^[4].re + a^[5].re;
  t2 := a^[4].im + a^[5].im;
  t3 := a^[6].re + a^[7].re;
  t4 := a^[6].im + a^[7].im;

  a^[5].re := a^[4].re - a^[5].re;
  a^[5].im := a^[4].im - a^[5].im;
  a^[7].re := a^[6].re - a^[7].re;
  a^[7].im := a^[6].im - a^[7].im;

  t5 := t2 - t4;
  t6 := t3 - t1;
  t1 := t1 + t3;
  t2 := t2 + t4;

  a^[6].re := a^[2].re - t5;
  a^[6].im := a^[2].im - t6;
  a^[4].re := a^[0].re - t1;
  a^[4].im := a^[0].im - t2;
  a^[2].re := a^[2].re + t5;
  a^[2].im := a^[2].im + t6;
  a^[0].re := a^[0].re + t1;
  a^[0].im := a^[0].im + t2;

  t1 := (a[5].re + a^[5].im) * SqrtHalf;
  t4 := (a[7].im + a^[7].re) * SqrtHalf;
  t2 := (a[5].im - a^[5].re) * SqrtHalf;
  t3 := (a[7].re - a^[7].im) * SqrtHalf;

  t6 := t3 - t1;
  t1 := t1 + t3;
  t5 := t2 - t4;
  t2 := t2 + t4;

  a^[7].re := a^[3].re - t5;
  a^[7].im := a^[3].im - t6;
  a^[5].re := a^[1].re - t1;
  a^[5].im := a^[1].im - t2;
  a^[3].re := a^[3].re + t5;
  a^[3].im := a^[3].im + t6;
  a^[1].re := a^[1].re + t1;
  a^[1].im := a^[1].im + t2;
end;

procedure u16(a: PComplexArray);
begin
  u8(a);
  u4(@a^[8]);
  u4(@a^[12]);

  UNTRANSFORMZERO(@a^[0], @a^[4], @a^[8], @a^[12]);
  UNTRANSFORMHALF(@a^[2], @a^[6], @a^[10], @a^[14]);
  UNTRANSFORM(@a^[1], @a^[5], @a^[9], @a^[13], d16[0].re, d16[0].im);
  UNTRANSFORM(@a^[3], @a^[7], @a^[11], @a^[15], d16[0].im, d16[0].re);
end;

procedure upass(a, w: PComplexArray; n: integer);
var
  k: integer;
  z: integer;
begin
  UNTRANSFORMZERO(@a^[0], @a^[2 * n], @a^[4 * n], @a^[6 * n]);
  k := 2 * n - 1;
  z := 0;

  while k > 0 do
  begin
    UNTRANSFORM(@a^[1 + z], @a^[2 * n + 1 + z], @a^[4 * n + 1 + z], @a^[6 * n + 1 + z], w^[z].re, w^[z].im);
    inc(z);
    dec(k);
  end;
end;

procedure u32(a: PComplexArray);
begin
  u16(a);
  u8(@a^[16]);
  u8(@a^[24]);
  upass(a, @d32, 4);
end;

procedure u64(a: PComplexArray);
begin
  u32(a);
  u16(@a^[32]);
  u16(@a^[48]);
  upass(a, @d64, 8);
end;

procedure u128(a: PComplexArray);
begin
  u64(a);
  u32(@a^[64]);
  u32(@a^[96]);
  upass(a, @d128, 16);
end;

procedure u256(a: PComplexArray);
begin
  u128(a);
  u64(@a^[128]);
  u64(@a^[192]);
  upass(a, @d256, 32);
end;

procedure u512(a: PComplexArray);
begin
  u256(a);
  u128(@a^[256]);
  u128(@a^[384]);
  upass(a, @d512, 64);
end;

procedure upassbig(a, w: PComplexArray; n: integer);
var
  k: integer;
  z: integer;
begin
  UNTRANSFORMZERO(@a^[0], @a^[2 * n], @a^[4 * n], @a^[6 * n]);
  k := n - 1;

  z := 0;
  while k > 0 do
  begin
    UNTRANSFORM(@a^[1 + z], @a^[2 * n + 1 + z], @a^[4 * n + 1 + z], @a^[6 * n + 1 + z], w^[z].re, w^[z].im);
    inc(z);
    dec(k);
  end;

  UNTRANSFORMHALF(@a^[1 + z], @a^[2 * n + 1 + z], @a^[4 * n + 1 + z], @a^[6 * n + 1 + z]);
  inc(z);
  k := n - 1;

  while k > 0 do
  begin
    UNTRANSFORM(@a^[1 + z], @a^[2 * n + 1 + z], @a^[4 * n + 1 + z], @a^[6 * n + 1 + z], w^[k - 1].im, w^[k - 1].re);
    inc(z);
    dec(k);
  end;
end;

procedure u1024(a: PComplexArray);
begin
  u512(a);
  u256(@a^[512]);
  u256(@a^[768]);
  upassbig(a, @d1024, 128);
end;

procedure u2048(a: PComplexArray);
begin
  u1024(a);
  u512(@a^[1024]);
  u512(@a^[1536]);
  upassbig(a, @d2048, 256);
end;

procedure u4096(a: PComplexArray);
begin
  u2048(a);
  u1024(@a^[2048]);
  u1024(@a^[3072]);
  upassbig(a, @d4096, 512);
end;

procedure u8192(a: PComplexArray);
begin
  u4096(a);
  u2048(@a^[4096]);
  u2048(@a^[6144]);
  upassbig(a, @d8192, 1024);
end;

procedure u16384(a: PComplexArray);
begin
  u8192(a);
  u4096(@a^[8192]);
  u4096(@a^[12288]);
  upassbig(a, @d16384, 2048);
end;

procedure fft2(a: PComplexArray);
begin
  c2(a);
end;

procedure fft4(a: PComplexArray);
begin
  c4(a);
end;

procedure fft8(a: PComplexArray);
begin
  c8(a);
end;

procedure fft16(a: PComplexArray);
begin
  c16(a);
end;

procedure fft32(a: PComplexArray);
begin
  c32(a);
end;

procedure fft64(a: PComplexArray);
begin
  c64(a);
end;

procedure fft128(a: PComplexArray);
begin
  c128(a);
end;

procedure fft256(a: PComplexArray);
begin
  c256(a);
end;

procedure fft512(a: PComplexArray);
begin
  c512(a);
end;

procedure fft1024(a: PComplexArray);
begin
  c1024(a);
end;

procedure fft2048(a: PComplexArray);
begin
  c2048(a);
end;

procedure fft4096(a: PComplexArray);
begin
  c4096(a);
end;

procedure fft8192(a: PComplexArray);
begin
  c8192(a);
end;

procedure fft16384(a: PComplexArray);
begin
  c16384(a);
end;

procedure ifft2(a: PComplexArray);
begin
  c2(a);
end;

procedure ifft4(a: PComplexArray);
begin
  u4(a);
end;

procedure ifft8(a: PComplexArray);
begin
  u8(a);
end;

procedure ifft16(a: PComplexArray);
begin
  u16(a);
end;

procedure ifft32(a: PComplexArray);
begin
  u32(a);
end;

procedure ifft64(a: PComplexArray);
begin
  u64(a);
end;

procedure ifft128(a: PComplexArray);
begin
  u128(a);
end;

procedure ifft256(a: PComplexArray);
begin
  u256(a);
end;

procedure ifft512(a: PComplexArray);
begin
  u512(a);
end;

procedure ifft1024(a: PComplexArray);
begin
  u1024(a);
end;

procedure ifft2048(a: PComplexArray);
begin
  u2048(a);
end;

procedure ifft4096(a: PComplexArray);
begin
  u4096(a);
end;

procedure ifft8192(a: PComplexArray);
begin
  u8192(a);
end;

procedure ifft16384(a: PComplexArray);
begin
  u16384(a);
end;

procedure scalec(a: PComplexArray; n: integer; u: Float);
var
  i: integer;
begin
  for i := 0 to n - 1 do
  begin
    a^[i].re := a^[i].re * u;
    a^[i].im := a^[i].im * u;
  end;
end;

procedure scale2(a: PComplexArray);
begin
  scalec(a, 2, scalec2);
end;

procedure scale4(a: PComplexArray);
begin
  scalec(a, 4, scalec4);
end;

procedure scale8(a: PComplexArray);
begin
  scalec(a, 8, scalec8);
end;

procedure scale16(a: PComplexArray);
begin
  scalec(a, 16, scalec16);
end;

procedure scale32(a: PComplexArray);
begin
  scalec(a, 32, scalec32);
end;

procedure scale64(a: PComplexArray);
begin
  scalec(a, 64, scalec64);
end;

procedure scale128(a: PComplexArray);
begin
  scalec(a, 128, scalec128);
end;

procedure scale256(a: PComplexArray);
begin
  scalec(a, 256, scalec256);
end;

procedure scale512(a: PComplexArray);
begin
  scalec(a, 512, scalec512);
end;

procedure scale1024(a: PComplexArray);
begin
  scalec(a, 1024, scalec1024);
end;

procedure scale2048(a: PComplexArray);
begin
  scalec(a, 2048, scalec2048);
end;

procedure scale4096(a: PComplexArray);
begin
  scalec(a, 4096, scalec4096);
end;

procedure scale8192(a: PComplexArray);
begin
  scalec(a, 8192, scalec8192);
end;

procedure scale16384(a: PComplexArray);
begin
  scalec(a, 16384, scalec16384);
end;

procedure mulc(a, b: PComplexArray; n: integer);
var
  t1, t2, t3, t4: Float;
  z: integer;
begin
  z := 0;
  while (z < n) do
  begin
    t1 := (a^[z].re * b^[z].re) - (a^[z].im * b^[z].im);
    t2 := (a^[z].im * b^[z].re) + (a^[z].re * b^[z].im);
    t3 := (a^[z + 1].re * b^[z + 1].re) - (a^[z + 1].im * b^[z + 1].im);
    t4 := (a^[z + 1].im * b^[z + 1].re) + (a^[z + 1].re * b^[z + 1].im);
    a^[z].re := t1;
    a^[z + 1].re := t3;
    a^[z].im := t2;
    a^[z + 1].im := t4;
    inc(z, 2);
  end;
end;

procedure mul2(a, b: PComplexArray);
begin
  mulc(a, b, 2);
end;

procedure mul4(a, b: PComplexArray);
begin
  mulc(a, b, 4);
end;

procedure mul8(a, b: PComplexArray);
begin
  mulc(a, b, 8);
end;

procedure mul16(a, b: PComplexArray);
begin
  mulc(a, b, 16);
end;

procedure mul32(a, b: PComplexArray);
begin
  mulc(a, b, 32);
end;

procedure mul64(a, b: PComplexArray);
begin
  mulc(a, b, 64);
end;

procedure mul128(a, b: PComplexArray);
begin
  mulc(a, b, 128);
end;

procedure mul256(a, b: PComplexArray);
begin
  mulc(a, b, 256);
end;

procedure mul512(a, b: PComplexArray);
begin
  mulc(a, b, 512);
end;

procedure mul1024(a, b: PComplexArray);
begin
  mulc(a, b, 1024);
end;

procedure mul2048(a, b: PComplexArray);
begin
  mulc(a, b, 2048);
end;

procedure mul4096(a, b: PComplexArray);
begin
  mulc(a, b, 4096);
end;

procedure mul8192(a, b: PComplexArray);
begin
  mulc(a, b, 8192);
end;

procedure mul16384(a, b: PComplexArray);
begin
  mulc(a, b, 16384);
end;

procedure InitSinCosTable;

  procedure SinCosTable(n: integer; s: PComplexArray);
  var
    tmp: Extended;
    i: integer;
    m: integer;
{$IFNDEF EXTENDED_PRECISION}
    fs, fc: Extended;
{$ENDIF}
  begin
    if n > 512 then
      m := n shr 3
    else
      m := n shr 2;
    tmp := 2 * PI / n;
    for i := 1 to m - 1 do
    begin
{$IFDEF EXTENDED_PRECISION}
      SinCos(tmp * i, s^[i - 1].im, s^[i - 1].re);
{$ELSE}
      SinCos(tmp * i, fs, fc);
      s^[i - 1].re := fc;
      s^[i - 1].im := fs;
{$ENDIF}
    end;
  end;

begin
  SinCosTable(16, @d16);
  SinCosTable(32, @d32);
  SinCosTable(64, @d64);
  SinCosTable(128, @d128);
  SinCosTable(256, @d256);
  SinCosTable(512, @d512);
  SinCosTable(1024, @d1024);
  SinCosTable(2048, @d2048);
  SinCosTable(4096, @d4096);
  SinCosTable(8192, @d8192);
  SinCosTable(16384, @d16384);
end;

procedure InitRevBitsTable;

  procedure CopyFromTo(pFrom, pTo: PWORDArray; x, z, m: PInteger);
  var
    i: integer;
  begin
    if x^ = 0 then
      x^ := 1
    else
      x^ := 0;
    z^ := z^ shl 1 - x^;
    m^ := m^ shl 1;
    for i := 0 to (m^ shr 1 - 1) do
    begin
      pTo^[i] := pFrom^[i] * 2;
      pTo^[i + m^ shr 1] := pFrom^[i] * 2;
    end;
    for i := z^ to (m^ - (m^ shr 1 - (z^ - 1))) do
      inc(pTo^[i]);
  end;

var
  i: integer;
  z: integer;
  x: integer;
  m: integer;
begin
  x := 1;
  z := 1;
  m := 2;
  for i := 0 to 1 do
    RevBits2[i] := i;
  CopyFromTo(@RevBits2, @RevBits4, @x, @z, @m);
  CopyFromTo(@RevBits4, @RevBits8, @x, @z, @m);
  CopyFromTo(@RevBits8, @RevBits16, @x, @z, @m);
  CopyFromTo(@RevBits16, @RevBits32, @x, @z, @m);
  CopyFromTo(@RevBits32, @RevBits64, @x, @z, @m);
  CopyFromTo(@RevBits64, @RevBits128, @x, @z, @m);
  CopyFromTo(@RevBits128, @RevBits256, @x, @z, @m);
  CopyFromTo(@RevBits256, @RevBits512, @x, @z, @m);
  CopyFromTo(@RevBits512, @RevBits1024, @x, @z, @m);
  CopyFromTo(@RevBits1024, @RevBits2048, @x, @z, @m);
  CopyFromTo(@RevBits2048, @RevBits4096, @x, @z, @m);
  CopyFromTo(@RevBits4096, @RevBits8192, @x, @z, @m);
  CopyFromTo(@RevBits8192, @RevBits16384, @x, @z, @m);
end;


initialization

InitRevBitsTable;
InitSinCosTable;

end.
