
$INSTALLDIR="c:\program files\Certara\nlme_engine"
if($INSTALLDIR -eq "" -or $INSTALLDIR -eq $null)
{
  "Installation directory is not specified"
  exit 1
}
powershell -noninteractive -executionpolicy remotesigned -File $INSTALLDIR\generic_run.ps1 none $INSTALLDIR $shared_directory D:\git\ADPOBenchmark\NLME\6\Run6 D:\git\ADPOBenchmark\NLME\6\Run6\jobControlFile.txt 1 WorkFlow
