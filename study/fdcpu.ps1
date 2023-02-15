#Powershell script for CPU Field Calculation

if ( $Args.Length -ne 17 ){
    Write-Host "Number of arguments is "$Args.Length" !!"
    Write-Host "We need 17 arguments"
    Write-Host "Usage: fdcpu Nx Ny Max_SNiter Max_ENiter Maxacc Gstartx Gbase Gth Overlap Diter EleFileName Plot start stop total RFV Numthread"
  exit 1
}

Write-Host "From $($Args[12]) to $($Args[13]) in total $($Args[14]) with RF $($Args[15]) V with $($Args[16]) Threads"
$i=$Args[12]
$V=0
$start  = date
while ( $i -lt $Args[13] ){
    $rad = $i / $Args[14] * 2 * [math]::pi
    $V= $Args[15] * [math]::cos( $rad )
    Write-Host "$i $V"
    mpiexec -n $Args[16] python E:\MyDocuments\Programs\cuda\pycudaseplap\LaplaceCylCorr.py $Args[0..9] $Args[11] $Args[10] $V
#    mpiexec -n $Args[16] python E:\MyDocuments\Programs\cuda\LaplaceCylElFuncFinal.py $Args[0..9] $Args[11] $Args[10] $V
    Copy-Item field.npy "$i-$($Args[14]).npy" -Force
    Remove-Item field.npy
    $i += 1
}
$end  = date
Write-Host $start $end
