$ErrorActionPreference = 'Stop'
if (-not $env:STATA_EXE) { throw 'Set STATA_EXE to your licensed Stata executable' }
$root = Split-Path $PSScriptRoot -Parent
$p = Start-Process -FilePath $env:STATA_EXE -ArgumentList '/e do tests/run.do' -WorkingDirectory $root -WindowStyle Hidden -PassThru
$p.WaitForExit()
if (-not (Select-String -Path (Join-Path $root 'run.log') -Pattern '^FOUR_PORT_TEST_EXECUTION_PASS$' -Quiet)) { throw 'Stata batch marker missing; inspect local log' }
python (Join-Path $PSScriptRoot 'compare.py')
if ($LASTEXITCODE -ne 0) { throw 'Frozen oracle comparison failed' }
