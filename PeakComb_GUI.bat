@echo off
setlocal EnableExtensions
chcp 65001 >nul

set "PROJECT_DIR=%~dp0"
set "PREFERRED_PYTHON=%USERPROFILE%\miniconda3\python.exe"

cd /d "%PROJECT_DIR%"

echo Starting PeakComb GUI...
echo Project: %CD%
echo.

if exist "%PREFERRED_PYTHON%" (
    echo Python: %PREFERRED_PYTHON%
    "%PREFERRED_PYTHON%" -m peakcomb
) else (
    where py >nul 2>nul
    if errorlevel 1 (
        python -m peakcomb
    ) else (
        py -3 -m peakcomb
    )
)

set "STATUS=%ERRORLEVEL%"
if not "%STATUS%"=="0" (
    echo.
    echo PeakComb exited with code %STATUS%.
    pause
)
exit /b %STATUS%
