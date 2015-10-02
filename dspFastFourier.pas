(* ********************************************************************
  *  dspFastFourier.pas                                               *
  *                                                                   *
  *  This unit is Part of the DC-DSP Component Pack v1.0              *
  *                                                                   *
  *  author    : Milenko Mitrovic                                     *
  *  email     : dcoder@dsp-worx.de                                   *
  *  web       : http://dsp-worx.de                                   *
  *  date      : 24-07-2003                                           *
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
  *  Based on D.J.Bernstein큦 Split Radix FFT C+ Code v0.76 (Linux).  *
  *  Available at http://cr.yp.to/djbfft.html.                        *
  *                                                                   *
  *  On my AthlonXP 1800+ CPU the FFT Results are the follwing :      *
  *  3DNowExt : FFT = 28% and IFFT = 24% faster then the FPU Code     *
  *  3DNow    : FFT = 18% and IFFT =  8% faster then the FPU Code     *
  *  SSE      : FFT = 17% and IFFT =  8% faster then the FPU Code     *
  *                                                                   *
  *  (C) 2003 Milenko Mitrovic <dcoder@dsp-worx.de>                   *
  *                                                                   *
  ******************************************************************** *)
{
  @abstract(Fast Fourier Transform Unit that can do Transformations
  up to a Size of 8192 Samples (<u>Power of 2 only!</u>).
  Optimized for SSE, 3DNow, 3DNowExt and CPUs that only supports FPU
  instructions. This Routine chooses automatically the fastest one.)
  @author(Milenko Mitrovic <dcoder@dsp-worx.de>)
  @created(Sep 01, 2002)
  @lastmod(Apr 02, 2004)
}

{$R-}
{$I dspFFTConfig.inc}
unit dspFastFourier;

interface

uses
  Classes, Math, dspConst, dspUtils, Windows, FastFourierDJB;

type
  { TDCFFT - Fast Fourier Transform Component that can do Transformations up to
    a Size of 8192 Samples (<u>Power of 2 only!</u>). Optimized for SSE, 3DNow,
    3DNowExt and CPUs that only supports FPU instructions. This Routine chooses
    automatically the fastest one. }
  TDCFFT = class(TComponent)
  private
    { @exclude }
    fSwapBuffer: Pointer;
    { @exclude }
    fCplx: PComplexArray;
    { @exclude }
    fScale: Boolean;
    { @exclude }
    fReOrder: Boolean;
    { @exclude }
    fFFTSize: TDCFFTSize;
    { @exclude }
    function GetFFTSize: integer;
  public
    { TDCFFT Constructor }
    constructor Create(AOwner: TComponent); override;
    { TDCFFT Destructor }
    destructor Destroy; override;
    { Does a forward Fourier Transform of the Samples in Complex. }
    procedure FFT;
    { Does a reverse Fourier Transform of the Samples in Complex. }
    procedure IFFT;
    { Pointer to an Array of 8192 Real and Imag Values that is used for the
      forward and inverse Transformation. This Buffer is Aligned on a 16 Byte
      boundary. }
    property Complex: PComplexArray read fCplx write fCplx;
    { Sets every Value to Zero. When creating the Component, this procedure is
      automatically called, so it큦 not needed to use it after creating. }
    procedure Flush;
  published
    { Normally this is needed for the Inverse Transform to Scaledown the Samples.
      Sometimes it큦 not needed. In that case set Scale to False. default = True }
    property Scale: Boolean read fScale write fScale;
    { Since this is a 2-2-4 Split Radix FFT, Samples are out of Order after the FFT.
      If Enabled then the Samples will be Reordered after the FFT. If ReOrder was
      Enabled for the forward FFT, then it MUST also be Enabled for the inverse FFT !
      Sometimes it큦 not needed to Reorder the Samples, so it can be disabled.
      default = True }
    property ReOrder: Boolean read fReOrder write fReOrder;
    { Sets the FFT Size that will be used for the Transformation. Can be one of
      TDCFFTSize Values. eg: fts2048 for a 2048 Point FFT. default = fts512 }
    property FFTSize: TDCFFTSize read fFFTSize write fFFTSize;
    { returns the Current FFT Size as an integer Value. }
    property FFTSizeInt: integer read GetFFTSize;
  end;

  { This procedure is made to use your own allocated Buffer to do a Fast Fourier
    Transform without using the TDCFFT Class. Values are the same as for the
    Component. When using a SSE capable CPU, make sure that your Buffer
    (Cplx:PComplexArray) is Aligned on 16 Byte Boundary or you will recieve an
    Error called "Privileged instruction". The Complex Buffer in the TDCFFT Class
    is already Aligned, so it큦 appricated to use this Class. SwapBuffer is used
    when ReOrder is used to prevent the temporary allocation/freeing of Memry used
    for the Reordering. Create a Buffer on your own, pass it to this function and
    Release it when the FFT is not needed anymore. If no Buffer is specified, the
    procedure will allocate the memory on it큦 own. For fastest Memory Movements
    this Buffer should be initialized using GetAlignedMemory. }
procedure dspDoFFT(Cplx: PComplexArray; FFTSize: integer; Inverse, Scale, ReOrder: Boolean; SwapBuffer: Pointer = nil);

implementation

uses SysUtils;

procedure scalec(a: PComplexArray; n: integer; u: Single);
asm
  push      ebx
  shl       edx, 3
  xor       ebx, ebx
  cmp       hasSSE, 0
  je        @@3DNOW
@@CSSE:
  cmp       edx, 16
  jl        @@3DNOW
  movss     xmm7, u
  shufps    xmm7, xmm7, $0
@@Loop2:
  movaps    xmm0, [eax+ebx]
  mulps     xmm0, xmm7
  movaps    [eax+ebx], xmm0
  add       ebx, 16
  cmp       ebx, edx
  jl        @@Loop2
  jmp       @@End
@@3DNOW:
  cmp       has3DNow, 0
  je        @@CFPU
  femms
  movd      mm7, u
  punpckldq mm7, mm7
@@Loop1:
  movq      mm0, [eax+ebx]
  pfmul     mm0, mm7
  movq      [eax+ebx], mm0
  add       ebx, 8
  cmp       ebx, edx
  jl        @@Loop1
  femms
  jmp       @@End
@@CFPU:
  fld       [eax+ebx]
  fmul      u
  fstp      [eax+ebx]
  add       ebx, 4
  cmp       ebx, edx
  jl        @@CFPU
@@End:
  pop       ebx
end;

procedure ReOrderFFT(a: PComplexArray; Size: integer; RevBin: PWORDArray; Inverse: Boolean;
  SwapBuffer: Pointer);
var
  tmp: PComplexArray;
  i: integer;
  FreeBuffer: Boolean;
begin
  if Assigned(SwapBuffer) then
  begin
    FreeBuffer := False;
    tmp := SwapBuffer;
  end
  else
  begin
    FreeBuffer := True;
    tmp := AllocMem(Size * SizeOf(TComplex));
  end;
  Move(a^, tmp^, Size * SizeOf(TComplex));
  if not Inverse then
    for i := 0 to Size - 1 do
      a[i] := tmp[RevBin[i]]
  else
    for i := 0 to Size - 1 do
      a[RevBin[i]] := tmp[i];
  if FreeBuffer then
    FreeMem(tmp);
end;

procedure dspDoFFT(Cplx: PComplexArray; FFTSize: integer; Inverse, Scale, ReOrder: Boolean; SwapBuffer: Pointer = nil);
begin
  if not Inverse then
  begin
    case FFTSize of
      2:
        begin
          fft2(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits2, Inverse, SwapBuffer);
        end;
      4:
        begin
          fft4(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits4, Inverse, SwapBuffer);
        end;
      8:
        begin
          fft8(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits8, Inverse, SwapBuffer);
        end;
      16:
        begin
          fft16(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits16, Inverse, SwapBuffer);
        end;
      32:
        begin
          fft32(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits32, Inverse, SwapBuffer);
        end;
      64:
        begin
          fft64(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits64, Inverse, SwapBuffer);
        end;
      128:
        begin
          fft128(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits128, Inverse, SwapBuffer);
        end;
      256:
        begin
          fft256(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits256, Inverse, SwapBuffer);
        end;
      512:
        begin
          fft512(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits512, Inverse, SwapBuffer);
        end;
      1024:
        begin
          fft1024(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits1024, Inverse, SwapBuffer);
        end;
      2048:
        begin
          fft2048(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits2048, Inverse, SwapBuffer);
        end;
      4096:
        begin
          fft4096(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits4096, Inverse, SwapBuffer);
        end;
      8192:
        begin
          fft8192(Cplx);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits8192, Inverse, SwapBuffer);
        end;
      16384:
        begin
          fft16384(Cplx);
          if ReOrder then ReOrderFFT(Cplx,FFTSize,@RevBits16384,Inverse,SwapBuffer);
        end;
    end;
  end
  else
  begin
    case FFTSize of
      2:
        begin
          if Scale then
            scalec(Cplx, 2, ScaleFFT2);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits2, Inverse, SwapBuffer);
          ifft2(Cplx);
        end;
      4:
        begin
          if Scale then
            scalec(Cplx, 4, ScaleFFT4);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits4, Inverse, SwapBuffer);
          ifft4(Cplx);
        end;
      8:
        begin
          if Scale then
            scalec(Cplx, 8, ScaleFFT8);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits8, Inverse, SwapBuffer);
          ifft8(Cplx);
        end;
      16:
        begin
          if Scale then
            scalec(Cplx, 16, ScaleFFT16);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits16, Inverse, SwapBuffer);
          ifft16(Cplx);
        end;
      32:
        begin
          if Scale then
            scalec(Cplx, 32, ScaleFFT32);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits32, Inverse, SwapBuffer);
          ifft32(Cplx);
        end;
      64:
        begin
          if Scale then
            scalec(Cplx, 64, ScaleFFT64);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits64, Inverse, SwapBuffer);
          ifft64(Cplx);
        end;
      128:
        begin
          if Scale then
            scalec(Cplx, 128, ScaleFFT128);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits128, Inverse, SwapBuffer);
          ifft128(Cplx);
        end;
      256:
        begin
          if Scale then
            scalec(Cplx, 256, ScaleFFT256);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits256, Inverse, SwapBuffer);
          ifft256(Cplx);
        end;
      512:
        begin
          if Scale then
            scalec(Cplx, 512, ScaleFFT512);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits512, Inverse, SwapBuffer);
          ifft512(Cplx);
        end;
      1024:
        begin
          if Scale then
            scalec(Cplx, 1024, ScaleFFT1024);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits1024, Inverse, SwapBuffer);
          ifft1024(Cplx);
        end;
      2048:
        begin
          if Scale then
            scalec(Cplx, 2048, ScaleFFT2048);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits2048, Inverse, SwapBuffer);
          ifft2048(Cplx);
        end;
      4096:
        begin
          if Scale then
            scalec(Cplx, 4096, ScaleFFT4096);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits4096, Inverse, SwapBuffer);
          ifft4096(Cplx);
        end;
      8192:
        begin
          if Scale then
            scalec(Cplx, 8192, ScaleFFT8192);
          if ReOrder then
            ReOrderFFT(Cplx, FFTSize, @RevBits8192, Inverse, SwapBuffer);
          ifft8192(Cplx);
        end;
      16384:
        begin
          if Scale   then scalec(Cplx,16384,ScaleFFT16384);
          if ReOrder then ReOrderFFT(Cplx,FFTSize,@RevBits16384,Inverse,SwapBuffer);
          ifft16384(Cplx);
        end;
    end;
  end;
end;

constructor TDCFFT.Create(AOwner: TComponent);
begin
  inherited Create(AOwner);

  fCplx := AllocMem(8192 * SizeOf(TComplex));
  fSwapBuffer := AllocMem(8192 * SizeOf(TComplex));
  Flush;
  fFFTSize := fts512;
  fReOrder := True;
  fScale := True;
end;

destructor TDCFFT.Destroy;
begin
  FreeMemory(fCplx);
  FreeMemory(fSwapBuffer);
  inherited Destroy;
end;

procedure TDCFFT.FFT;
begin
  dspDoFFT(fCplx, 1 shl (integer(fFFTSize) + 1), False, fScale, fReOrder, fSwapBuffer);
end;

procedure TDCFFT.IFFT;
begin
  dspDoFFT(fCplx, 1 shl (integer(fFFTSize) + 1), True, fScale, fReOrder, fSwapBuffer);
end;

function TDCFFT.GetFFTSize: integer;
begin
  Result := 1 shl (integer(fFFTSize) + 1);
end;

procedure TDCFFT.Flush;
begin
  FillChar(fCplx^, 8192 * SizeOf(TComplex), 0);
end;

end.
