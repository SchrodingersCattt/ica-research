@echo off
echo ============================================
echo  NCI/IRI Render Pipeline v2
echo  Chalk Matte + Wireframe Isosurface
echo ============================================
echo.

if "%ICA_VMD_EXE%"=="" set ICA_VMD_EXE=vmd.exe
if "%ICA_TACHYON_EXE%"=="" set ICA_TACHYON_EXE=tachyon_WIN32.exe
if "%ICA_NCI_ROOT%"=="" (
  echo Set ICA_NCI_ROOT to the local NCI input directory first.
  exit /b 2
)
set VMD=%ICA_VMD_EXE%
set TACHYON=%ICA_TACHYON_EXE%
set BASEDIR=%ICA_NCI_ROOT%
set TCL=%BASEDIR%\render_nci_v2_auto.tcl
set SCENE=%BASEDIR%\scene_v2.dat
set TGA=%BASEDIR%\NCI_3d_v2_quick.tga
set PNG=%BASEDIR%\NCI_3d_v2_quick.png

echo [1/3] Launching VMD to export scene...
%VMD% -e "%TCL%"
echo VMD done.
echo.

if not exist "%SCENE%" (
    echo ERROR: Scene file not found: %SCENE%
    echo VMD may have failed. Check for errors.
    pause
    exit /b 1
)

echo [2/3] Running Tachyon ray tracer (1200x1600, 4 AA samples)...
%TACHYON% "%SCENE%" -aasamples 4 -res 1200 1600 -format TGA -o "%TGA%"
echo Tachyon done.
echo.

if not exist "%TGA%" (
    echo ERROR: TGA file not found: %TGA%
    pause
    exit /b 1
)

echo [3/3] Converting TGA to PNG...
magick "%TGA%" "%PNG%"
echo.

if exist "%PNG%" (
    echo ============================================
    echo  SUCCESS! Output:
    echo    %PNG%
    echo ============================================
) else (
    echo PNG conversion failed. TGA available at:
    echo    %TGA%
    echo Try: magick "%TGA%" "%PNG%"
)

echo.
pause
