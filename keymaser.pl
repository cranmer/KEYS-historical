#!/usr/local/bin/perl 
#KEYS - Kernel Estimation of Your Shapes
#
#I am the keymaster, are you the gatekeeper?
#
#Author - Kyle Cranmer
#Date - July 1, 1999 
#
#Arguments - 	First argument is the inputfile 
#
#Documentation at: 
#    http://www-wisconsin.cern.ch/~cranmer/keys.html
#
$gatekeeper = shift;
($logfile, $bunk) = split /\./, $gatekeeper;
$logfile = "$logfile".".log";


#	##############################################
#		Task 1 - Input Control File
#	##############################################
open KEEPER, $gatekeeper or die "Can't open $gatekeeper\n";
open LOGFILE, ">$logfile";

$funcno = 0;
$fileno = 0;
$lookforfunction = 1;
$lookforfile = 1;
$cuts2 = 1;
chomp($line = <KEEPER>);
while ($cuts2) {
    if($line =~ /CUT/) {
	push @cuts, $line;
	chomp($line = <KEEPER>);
    }
    else {
	$cuts2 = 0;
    }
}

while($lookforfunction) {
    $funcno++;
    $fileno = 0;
    $lookforfile = 1;
#    print "Got to look for function $funcno, $fileno\n";
    ($bunk, $fu) = split /:\s+/, $line;
    ($fun,$bunk) = split /\s+/, $fu;
    push @func, $fun;
    ($funn, $bunk) = split /\./, $fun;
    push @funcname, $funn;
    chomp($line = <KEEPER>);
    ($bunk, $lu) = split /:\s+/, $line;
    ($lum, $bunk) = split /\s+/, $lu;
    push @lumi, $lum;
    chomp($line = <KEEPER>);
    ($bunk, $lb) = split /:\s+/, $line;
    ($lbn, $bunk) = split/\s+/, $lb;
    push @lbnd, $lbn;
    chomp($line = <KEEPER>);
    ($bunk, $sl) = split /:\s+/, $line;
    ($slb, $bunk) = split /\s+/, $sl;
    push @slbd, $slb;
    chomp($line = <KEEPER>);
    ($bunk, $ub) = split /:\s+/, $line;
    ($ubn, $bunk) = split /\s+/, $ub;
    push @ubnd, $ubn;
    chomp($line = <KEEPER>);
    ($bunk, $su) = split /:\s+/, $line;
    ($sub, $bunk) = split /\s+/, $su;
    push @subd, $sub;
    chomp($line = <KEEPER>);
    ($bunk, $so) = split /:\s+/, $line;
    ($smo, $bunk) = split /\s+/, $so;
    push @smoo, $smo;
    chomp($line = <KEEPER>);
    while($lookforfile){
	$fileno++;
	$nfiles[$funcno - 1] = $fileno;
#	print "Got to look for file $funcno, $fileno\n";
	($bunk, $nt) = split /:\s+/, $line;
	($ntu, $bunk) = split /\s+/, $nt;
	chomp($line = <KEEPER>);
	($bunk, $ni) = split /:\s+/, $line;
	($nti, $bunk) = split /\s+/, $ni;
	chomp($line = <KEEPER>);
	($bunk, $va) = split /:\s+/, $line;
	($var, $bunk) = split /\s+/, $va;
	chomp($line = <KEEPER>);
	($bunk, $se) = split /:\s+/, $line;
	($sel, $bunk) = split /\s+/, $se;
	chomp($line = <KEEPER>);
	($bunk, $no) = split /:\s+/, $line;
	($num, $bunk) = split /\s+/, $no;
	chomp($line = <KEEPER>);
	($bunk, $xs) = split /:\s+/, $line;
	($xse, $bunk) = split /\s+/, $xs;
	chomp($line = <KEEPER>);
	if ($line =~ /\s+Filename/) {$lookforfile = 1;}
	else {$lookforfile = 0;}
#	print "lookforfile = $lookforfile\n";
	$ntup[$funcno-1][$fileno-1] = $ntu;
	$ntid[$funcno-1][$fileno-1] = $nti;
	$vari[$funcno-1][$fileno-1] = $var;
	$self[$funcno-1][$fileno-1] = $sel;
	$numb[$funcno-1][$fileno-1] = $num;
	$norm[$funcno-1][$fileno-1] = $xse/$num;
	$xsec[$funcno-1][$fileno-1] = $xse;

    }
    if ($line =~ /OUTPUTFILE/) {$lookforfunction = 1;}
    else {$lookforfunction = 0;}
#    print "lookforfunction = $lookforfunction\n";
#    print "function $funcno\n";
#    print "nfiles = $nfiles[$funcno-1]\n";
    $NormSum = 0;
    for ($fi = 0; $fi < $nfiles[$funcno-1]; $fi++){
	$NormSum += $norm[$funcno-1][$fi];
	$NORMSUM[$funcno-1] = $NormSum;
    }
    for ($fi = 0; $fi < $nfiles[$funcno-1]; $fi++){
#	print "NormSum = $NormSum\n";
#	print LOGFILE "norm is $norm[$funcno-1][$fi]\n";
      	$weights[$funcno-1][$fi] = $norm[$funcno-1][$fi]/$NormSum;
#	print "Weights $weights[$funcno-1][$fi]\n";
    }
}
#print "$func[0]\n";
#print "$ntup[0][0]\n";
#print "$xsec[0][0]\n";
#print "$func[1]\n";
#print "$ntup[1][0]\n";
#print "$xsec[1][0]\n";
#print "$func[1]\n";
#print "$ntup[1][1]\n";
#print "$xsec[1][1]\n";


close KEEPER;

#	##############################################
#		Task 2 - Make Individual functions
#	##############################################

#For each Output File - run steps 1 - 5
for ($fn = 0; $fn < scalar(@func); $fn++) {

    print LOGFILE "Working on function $func[$fn]\n";
######################################################
#Step 1 - Make kumac to dump ntuple data, then invoke paw
    for ($fi = 0; $fi < $nfiles[$fn]; $fi++) {
#	print LOGFILE "Working on file $ntup[$fn][$fi]\n";
	$dump = "dump$fn"."_$fi".".kumac"; 
#	print $dump;
	open DUMP, ">$dump"  
	    or die "Can't open temporary kumac to get events from ntuple\n";
	foreach $thiscut(@cuts) {
#	    print DUMP $thiscut."\n";
	    ($bunk, $cutnumber, $cutexpression) = split /\s/, $thiscut;
	    print DUMP "CUT $cutnumber ($cutexpression)\n";
	}
	print DUMP "hi/file 1 $ntup[$fn][$fi]\n";
	print DUMP "nt/cut 98 $self[$fn][$fi]\n";
	print DUMP "nt/cut 97 $lbnd[$fn]".".le.$vari[$fn][$fi]\n";
	print DUMP "nt/cut 96 $ubnd[$fn]".".ge.$vari[$fn][$fi]\n";
	if($slbd[$fn] ==1 and $subd[$fn] ==1) {
	    print DUMP "nt/cut 99 \$96.and.\$97.and.\$98\n";
	}
	if($slbd[$fn] ==1 and $subd[$fn] !=1) {
	    print DUMP "nt/cut 99 \$97.and.\$98\n";
	}
	if($slbd[$fn] !=1 and $subd[$fn] ==1) {
	    print DUMP "nt/cut 99 \$96.and.\$98\n";
	}
	if($slbd[$fn] !=1 and $subd[$fn] !=1) {
	    print DUMP "nt/cut 99 \$98\n";
	}
	print DUMP "nt/dump $ntid[$fn][$fi]".
	    ".$vari[$fn][$fi] \$99 ! 1 'data$fn"."_$fi".".txt'\n";
       	print DUMP "close 1\n";
	close DUMP;
	($kumacname, @bunk2) = split /\./, $dump;
#	print "\n\n$kumacname\n\n";
	`paw -b $kumacname > delme.txt`;
	unlink "$dump" or warn "Can't clean up temporary kumacs\n";
    }


######################################################
#Step 2 - Input each ntuple dump file seperately
    $neventstot = 0;
    for ($fi = 0; $fi < $nfiles[$fn]; $fi++) {
	$eventsinfile = "data$fn"."_$fi".".txt";
#	print "length of data array = ".scalar(@data)."\n";;
#	print "inputfile = ".$eventsinfile."\n";
	open EVENTSIN, "$eventsinfile" or die "Can't open $eventsinfile\n";
	$nevents = 0;
	while(<EVENTSIN>){
	    $notempty = 1;
	    chomp($line = $_);
	    if($line =~ /\S/) {
		($number, @bunk) = split /\n/, $line;
		push @data, $number;
		$nevents++;
	    }
	}
# next if is to get last line, input is kinda wierd
#	if($line =~ /\S/) {
#	    ($number, @bunk) = split /\n/, $line;
#	    push @data, $number;
#	    $nevents++;
#	}

#	print "length of data array = ".scalar(@data)."\n";;
#	print "firstdatapoint is $data[0]\n";
	close INFILE;
	unlink $eventsinfile or warn "Can't clean up temporary data files\n";
#	$nevents--;
	if ($nevents > 0) {$neventstot += $nevents;}
	else {$nevents = 0;}
	$neventsinfile[$fn][$fi] = $nevents;
	print LOGFILE"Events in File $fi = $nevents\n";
    }
    if ($neventstot == 0) {
	$bins[$fn] = 1;
    }
    else {
	$bins[$fn] = 4*(1 + log($neventstot)/log(2.));
    }
    use Env qw(QBIN);
    if ($QBIN > 0) {
	$bins[$fn] = $QBIN;
    }
    else {
	$bins[$fn] = int $bins[$fn];
    }
    print LOGFILE "nbins =$bins[$fn]\n";
    print LOGFILE "Total Number of Events =$neventstot\n";
    $normforhist =($ubnd[$fn]-$lbnd[$fn])/$bins[$fn]*$lumi[$fn]*$NORMSUM[$fn];
    $binwidth = ($ubnd[$fn]-$lbnd[$fn])/$bins[$fn];
    $NPAReff = 0;
    for ($fii = 0; $fii < $nfiles[$fn]; $fii++) {
	$NPAReff += $neventsinfile[$fn][$fii]*$weights[$fn][$fii];
    }
    $normforhist *= $NPAReff;

######################################################
#Step 3 - Write out keys function 

#####################################################
#Part 1 - Make a function to output data points in increasing order
# So that the Kolmogorov-Smirnov test can be calculated.
#This is actually the input to another function b/c KS.F cannot have more than
#100 continued lines (about 300 data points)

    $shapeinputfile = "$funcname[$fn].dat";
    open SHAPEINPUT, ">$shapeinputfile" or die "Can't open perl created input data file\n";
    foreach $dataevent (@data) {
	print SHAPEINPUT "$dataevent\n";
    }
    close SHAPEINPUT;

    foreach $event  (@data) {
	($value, $exponent) = split /e/, $event;
	push @temp, $value*10**$exponent;
    }
    sub numerically {$a <=> $b;}
    @sorted_data = sort numerically @temp;

    $ksshapeinputfile = "ks$funcname[$fn].dat";
    open KSSHAPEINPUT, ">$ksshapeinputfile" or die "Can't open perl created input data file\n";
    foreach $dataevent (@sorted_data) {
	print KSSHAPEINPUT "$dataevent\n";
    }
    close KSSHAPEINPUT;

#####################################################
#Part 2 - Make a function to evaluate the true Keys function 
# The true Keys function is too slow, so we evaluate it at 300 pts
# then use that to plot it in Paw and to generate a Linear Interpolated 
# approximation which is the final product

    $shapefile = "shape.F";
    open SHAPE, ">$shapefile" or die "Can't create file for plotting\n";
    print SHAPE "      PROGRAM shape\n";
    print SHAPE "      EXTERNAL $funcname[$fn]\n";
    print SHAPE "      REAL $funcname[$fn]\n";
    print SHAPE "      REAL X,Y, NORMFORHIST, LBOUND, UBOUND, RANGE\n";
    $range = $ubnd[$fn] - $lbnd[$fn];
    print SHAPE "      DATA LBOUND, UBOUND, RANGE /$lbnd[$fn] , $ubnd[$fn], $range /\n";
    print SHAPE "      DO 200 Q = 0, 300\n";
        print SHAPE "         X = $range"."\*Q/300. + $lbnd[$fn]\n";
    print SHAPE "         Y = $funcname[$fn]"."(X)\*$normforhist \n";
    print SHAPE "         WRITE(\*,\*) Y\n";
    print SHAPE "  200 CONTINUE\n";
    print SHAPE "      END\n";
    close SHAPE;

#####################################################
#Part 3 - Make a function to evaluate the KS.F function for input into Paw 
# This is necessary because it avoids problems with large Data Sets.
# Paw cannot interpret fortran that reads from files and it there is a limit 
# on the number of continued lines.

    $ksshapefile = "ksshape.F";
    open KSSHAPE, ">$ksshapefile" or die "Can't create file for plotting\n";
    print KSSHAPE "      PROGRAM KSshape\n";
    print KSSHAPE "      EXTERNAL ks\n";
    print KSSHAPE "      REAL X,Y, KS\n";
    print KSSHAPE "      DO 200 Q = 0, 300\n";
    print KSSHAPE "         X = $range"."\*Q/300. + $lbnd[$fn]\n";
    print KSSHAPE "         Y = ks"."(X)\*1.\n";
    print KSSHAPE "         WRITE(\*,\*) Y\n";
    print KSSHAPE "  200 CONTINUE\n";
    print KSSHAPE "      END\n";
    close KSSHAPE;


#####################################################
#Part 4 - Make the real True Keys function
# This function gets compiled, evaluated, and approximated to produce the 
# final function.  It reads the events from an input file b/c the limit on 
# the number of continued lines.  This function is pretty slow - it has to 
# evaluate N^2 Gaussians.  Some signal shapes have ~8000 events, that's 
# 64,000,000 Gaussians per point that's evaluated!


$outfile = "keys$func[$fn]";
open OUTF, ">$outfile" or die "Can't create $func[$fn]\n";

print OUTF "      REAL FUNCTION $funcname[$fn]"."(V)\n";
print OUTF "      EXTERNAL g\n";
print OUTF "      REAL V( 1)\n";
print OUTF "      REAL g\n";
print OUTF "      common /cb_$funcname[$fn] /\n";
print OUTF "     +h, ISTHISPI, SIGMA2, EVENTS, WEIGHTS, \n";
print OUTF "     +LBOUND, UBOUND, NPAR,hs\n";
print OUTF "      REAL h\n";
print OUTF "      REAL ISTHISPI\n";
print OUTF "      REAL SIGMA2\n";
print OUTF "      INTEGER NPAR\n";
print OUTF "      REAL EVENTS($neventstot)\n";
print OUTF "      REAL WEIGHTS($neventstot)\n";
print OUTF "      REAL LBOUND, UBOUND\n";
print OUTF "      INTEGER NFILES\n";
print OUTF "      REAL f\n";
print OUTF "      REAL NORM($nfiles[$fn])\n";
print OUTF "      REAL NORMSUM\n";
print OUTF "      INTEGER FIRSTWITHNORM($nfiles[$fn])\n";
print OUTF "      LOGICAL FIRST\n";
print OUTF "      SAVE FIRST\n";
print OUTF "      DATA  FIRST / .true. /\n";
print OUTF "      REAL hs($neventstot), a, y, yi, Ni, mirrorU, mirrorL, Alpha\n";
print OUTF "      REAL MEAN, hsigma\n";
print OUTF "      REAL NORMFORHIST\n";
print OUTF "      DATA NORMFORHIST, Alpha / $normforhist , $smoo[$fn] /\n";
print OUTF "      DATA NPAR, NFILES/ $neventstot, $nfiles[$fn] /\n";
print OUTF "      DATA MEAN, SIGMA2/ 0., 0. /\n";
print OUTF "      DATA LBOUND, UBOUND, ISTHISPI / $lbnd[$fn], $ubnd[$fn], 3.1415926535 /\n";


#this is the code to put the data in the function itself instead of a read sta.
#replaced with read statement from file
#print OUTF "      DATA EVENTS / ";
#$j=0;
#for ($j=0; $j<$neventstot; $j++){
#	if($j==0){
#		print OUTF "$data[$j]";
#		if($j<($neventstot-1)){
#		    $j++;
#		    print OUTF ", $data[$j]";
#		}
#		print OUTF "\n";
#	}
#	else{
#		print OUTF "     +, $data[$j]";
#		if($j<($neventstot-1)){
#		    $j++;
#		    print OUTF ", $data[$j]";
#		}
#		if($j<($neventstot-1)){
#		    $j++;
#		    print OUTF ", $data[$j]";
#		}
#		print OUTF "\n";
#	}
#} 
#print OUTF "     +/\n";

print OUTF "      DATA NORM / $norm[$fn][0]\n";
for ($fi=1; $fi<$nfiles[$fn] ; $fi++) {
    print OUTF "     +, $norm[$fn][$fi]\n";
}
print OUTF "     +/\n";
print OUTF "      DATA FIRSTWITHNORM / 1";
    $first = 1;
for ($fi=1; $fi<$nfiles[$fn] ; $fi++) {
    $first += $neventsinfile[$fn][$fi - 1];
    print OUTF ", $first ";
}
print OUTF "/\n";

print OUTF "      OPEN(unit=22, file = \n";
print OUTF "     +'"."$shapeinputfile"."'\n";
print OUTF "     +, status = 'old')\n";
print OUTF "      DO 22 i = 1, NPAR\n";
print OUTF "      read(22,\*), EVENTS(i)\n";
print OUTF "  22  CONTINUE\n";
print OUTF "      close(22)\n";

print OUTF "      h = ((4/3)**(.2))*(NPAR**(-.2))\n";
print OUTF "      NORMSUM = 0.\n";
print OUTF "      DO 10 HH = 1, (NFILES-1)\n";
print OUTF "      NORMSUM = NORMSUM + \n";
print OUTF "     +(FIRSTWITHNORM(HH+1)-FIRSTWITHNORM(HH))*NORM(HH)\n";
print OUTF "   10 CONTINUE\n";
print OUTF "      NORMSUM = NORMSUM + \n";
print OUTF "     +(NPAR-FIRSTWITHNORM(NFILES)+1)*NORM(NFILES)\n";

print OUTF "      DO 11 FF = 1, (NFILES-1)\n";
print OUTF "      DO 12 GG = FIRSTWITHNORM(FF), (FIRSTWITHNORM(FF+1)-1)\n";
print OUTF "      WEIGHTS(GG) = NORM(FF)/NORMSUM\n";
print OUTF "   12 CONTINUE\n";
print OUTF "   11 CONTINUE\n";

print OUTF "      DO 13 D = FIRSTWITHNORM(NFILES), NPAR\n";
print OUTF "      WEIGHTS(D) = NORM(NFILES)/NORMSUM\n";
print OUTF "   13 CONTINUE\n";

print OUTF "      MEAN = 0.\n"; 
print OUTF "      DO 15 I = 1, NPAR\n";
print OUTF "      MEAN = MEAN + WEIGHTS(I)*EVENTS(I)\n";
print OUTF "   15 CONTINUE\n";

print OUTF "      SIGMA2 = 0.\n";
print OUTF "      DO 20 J = 1, NPAR\n";
print OUTF "      SIGMA2 = SIGMA2 + WEIGHTS(J)*(MEAN - EVENTS(J))**2\n";
print OUTF "   20 CONTINUE\n";
print OUTF "      f = 0.\n";

    if ($slbd[$fn]==1) {
print OUTF "      IF (V(1).lt.LBOUND) THEN\n";
print OUTF "      DO 30 L = 1, NPAR\n";
print OUTF "      IF (EVENTS(L).lt.LBOUND) THEN\n";
print OUTF "      f = f + 0\n";
print OUTF "      END IF\n";
print OUTF "   30 CONTINUE\n";
print OUTF "      ELSE";
}
    if ($subd[$fn]==1) {
print OUTF "      IF (V(1).gt.UBOUND) THEN\n";
print OUTF "      DO 40 M = 1, NPAR\n";
print OUTF "      IF (EVENTS(M).gt.UBOUND) THEN\n";
print OUTF "      f = f + 0\n";
print OUTF "      END IF\n";
print OUTF "   40 CONTINUE\n";
print OUTF "      ELSE\n";
}
print OUTF "\n";
print OUTF "      DO 50 K = 1, NPAR\n";
print OUTF "      if(first)then\n";
print OUTF "      hs(k) = (h*Sqrt(Sqrt(Sigma2)/g(EVENTS(K))))*\n";
print OUTF "     + Alpha/(2.0*Sqrt(3.0))\n";
print OUTF "      hsigma = h*Sqrt(Sigma2)\n"; 
print OUTF "      If(hs(k).lt.(hsigma/10.0)) THEN\n";
print OUTF "      hs(k) = hsigma/10.0\n";
print OUTF "      END IF\n";
print OUTF "      END IF\n";
#print OUTF "*      hs = h*Sqrt(Sigma2)\n";
#    if($slbd[$fn]==1 and $subd[$fn]==1){
#print OUTF "      IF ((EVENTS(K).eq.LBOUND).or.(EVENTS(K).eq.UBOUND)) THEN\n";
#print OUTF "      GOTO 50\n";
#print OUTF "      ELSE";
#}
#    elsif($slbd[$fn]==1){
#print OUTF "      IF (EVENTS(K).eq.LBOUND) THEN\n";
#print OUTF "      GOTO 50\n";
#print OUTF "      ELSE";
#    }
#    elsif($subd[$fn]==1){
#print OUTF "      IF (EVENTS(K).eq.UBOUND) THEN\n";
#print OUTF "      GOTO 50\n";
#print OUTF "      ELSE";
#    }
    if ($slbd[$fn]==1){
print OUTF "      IF (EVENTS(K).lt.(LBOUND+2*hs(k))) THEN\n";
print OUTF "      mirrorL = 2*LBOUND-Events(K)\n";
print OUTF "      f = f + WEIGHTS(K)*EXP((-(V(1)-MirrorL)**2)/(2.0*hs(k)*hs(k)))\n";
print OUTF "     +/(hs(k)*SQRT(2.0*IsThisPi))\n";
print OUTF "      f = f +WEIGHTS(K)*EXP((-(V(1)-EVENTS(K))**2)/(2.0*hs(k)*hs(k)))\n";
print OUTF "     +/(hs(k)*SQRT(2.0*IsThisPi))\n";

print OUTF "      ELSE";
}
    if($subd[$fn]==1){
print OUTF "      IF (EVENTS(K).gt.(UBOUND-2*hs(k))) THEN\n";
print OUTF "      mirrorU = 2*UBOUND-Events(K)\n";
print OUTF "      f = f + WEIGHTS(K)*EXP((-(V(1)-MirrorU)**2)/(2.0*hs(k)*hs(k)))\n";
print OUTF "     +/(hs(k)*SQRT(2.0*IsThisPi))\n";
print OUTF "      f = f +WEIGHTS(K)*EXP((-(V(1)-EVENTS(K))**2)/(2.0*hs(k)*hs(k)))\n";
print OUTF "     +/(hs(k)*SQRT(2.0*IsThisPi))\n";

print OUTF "      ELSE\n";
}
print OUTF "\n";
print OUTF "      f = f +WEIGHTS(K)*EXP((-(V(1)-EVENTS(K))**2)/(2.0*hs(k)*hs(k)))\n";
print OUTF "     +/(hs(k)*SQRT(2.0*IsThisPi))\n";

    if($slbd[$fn]==1 or $subd[$fn]==1) {
print OUTF "      END IF\n";
}
print OUTF "   50 CONTINUE\n";
    if($slbd[$fn]==1 or $subd[$fn]==1) {
print OUTF "      END IF\n";
}
print OUTF "      $funcname[$fn] = f\n";
print OUTF "      first= .false.\n";
print OUTF "      RETURN\n";
#print OUTF "      f = f\n";
#print OUTF "      f = NORMFORHIST*f\n";
print OUTF "      END\n";



print OUTF "      REAL FUNCTION g(cV)\n";
print OUTF "      common /cb_$funcname[$fn] /\n";
print OUTF "     +h, IsThisPi, SIGMA2, EVENTS, WEIGHTS, \n";
print OUTF "     +LBOUND, UBOUND, NPAR,dum\n";
print OUTF "      REAL h\n";
print OUTF "      REAL ISTHISPI\n";
print OUTF "      REAL SIGMA2\n";
print OUTF "      INTEGER NPAR\n";
print OUTF "      REAL EVENTS($neventstot)\n";
print OUTF "      REAL dum($neventstot)\n";
print OUTF "      REAL WEIGHTS($neventstot)\n";
print OUTF "      REAL LBOUND, UBOUND\n";
print OUTF "      REAL MirrorL, MirrorU\n";
print OUTF "      REAL cV\n";
print OUTF "      g = 0.\n";
print OUTF "      hs = h*SQRT(SIGMA2)\n";

    if($slbd[$fn]==1){
print OUTF "      IF (cV.lt.LBOUND) THEN\n";
print OUTF "      DO 30 L = 1, NPAR\n";
print OUTF "      IF (EVENTS(L).eq.LBOUND) THEN\n";
print OUTF "      g = g + 1\n";
print OUTF "      END IF\n";
print OUTF "   30 CONTINUE\n";
print OUTF "      ELSE";
}
    if($subd[$fn]==1){
print OUTF "      IF (cV.gt.UBOUND) THEN\n";
print OUTF "      DO 40 M = 1, NPAR\n";
print OUTF "      IF (EVENTS(M).eq.UBOUND) THEN\n";
print OUTF "      g = g + 1\n";
print OUTF "      END IF\n";
print OUTF "   40 CONTINUE\n";
print OUTF "      ELSE\n";
}
print OUTF "\n";
print OUTF "      DO 50 K = 1, NPAR\n";
#    if($slbd[$fn]==1 and $subd[$fn]==1){
#print OUTF "      IF ((EVENTS(K).eq.LBOUND).or.(EVENTS(K).eq.UBOUND)) THEN\n";
#print OUTF "      GOTO 50\n";
#print OUTF "      ELSE";
#}
#    elsif($slbd[$fn]==1){
#print OUTF "      IF (EVENTS(K).eq.LBOUND) THEN\n";
#print OUTF "      GOTO 50\n";
#print OUTF "      ELSE";
#    }
#    elsif($subd[$fn]==1){
#print OUTF "      IF (EVENTS(K).eq.UBOUND) THEN\n";
#print OUTF "      GOTO 50\n";
#print OUTF "      ELSE";
#    }
    if($slbd[$fn]==1){
print OUTF "      IF (EVENTS(K).lt.(LBOUND+2*hs)) THEN\n";
print OUTF "      mirrorL = 2*LBOUND-Events(K)\n";
print OUTF "      g = g + WEIGHTS(K)*EXP((-(cV-MirrorL)**2)/(2.0*hs*hs))\n";
print OUTF "     +/(hs*SQRT(2.0*IsThisPi))\n";
print OUTF "      g = g + WEIGHTS(K)*EXP((-(cV-EVENTS(K))**2)/(2.0*hs*hs))\n";
print OUTF "     +/(hs*SQRT(2.0*IsThisPi))\n";

print OUTF "      ELSE";
}
    if($subd[$fn]==1){
print OUTF "      IF (EVENTS(K).gt.(UBOUND-2*hs)) THEN\n";
print OUTF "      mirrorU = 2*UBOUND-Events(K)\n";
print OUTF "      g = g + WEIGHTS(K)*EXP((-(cV-MirrorU)**2)/(2.0*hs*hs))\n";
print OUTF "     +/(hs*SQRT(2.0*IsThisPi))\n";
print OUTF "      g = g + WEIGHTS(K)*EXP((-(cV-EVENTS(K))**2)/(2.0*hs*hs))\n";
print OUTF "     +/(hs*SQRT(2.0*IsThisPi))\n";

print OUTF "      ELSE\n";
}
print OUTF "\n";
print OUTF "      g = g + WEIGHTS(K)*EXP((-(cV-EVENTS(K))**2)/(2.0*hs*hs))\n";
print OUTF "     +/(hs*SQRT(2.0*IsThisPi))\n";
    if($slbd[$fn]==1 or $subd[$fn]==1) {
print OUTF "      END IF\n";
}
print OUTF "   50 CONTINUE\n";
    if($slbd[$fn]==1 or $subd[$fn]==1) {
print OUTF "      END IF\n";
}
print OUTF "      g = g\n";
print OUTF "      RETURN\n";
print OUTF "      END\n";

close OUTF;

######################################################
#Step 4 - Write out Data's cumulative distribution function 


$ksfor = "ks.F";
open KS, ">$ksfor" or die "Can't Open $ksfor";
print KS "      REAL FUNCTION ks(V)\n";
print KS "      REAL V(1)\n";
print KS "      REAL EVENTS($neventstot)\n";
print KS "      INTEGER NPAR\n";
print KS "      DATA NPAR / $neventstot /\n";

print KS "      OPEN(unit=22, file = \n";
print KS "     +'"."$ksshapeinputfile"."'\n";
print KS "     +, status = 'old')\n";
print KS "      DO 222 ii = 1, NPAR\n";
print KS "      read(22,\*), EVENTS(ii)\n";
print KS " 222  CONTINUE\n";
print KS "      close(22)\n";


#    print KS "      DATA EVENTS / ";
#    $j=0;
#    for ($j=0; $j<$neventstot; $j++){
#	if($j==0){
#	    print KS "$sorted_data[$j]";
#	    if($j<($neventstot-1)){
#		$j++;
#		print KS ", $sorted_data[$j]";
#	    }
#	    print KS "\n";
#	}
#	else{
#	    print KS "     +, $sorted_data[$j]";
#	    if($j<($neventstot-1)){
#		$j++;
#		print KS ", $sorted_data[$j]";
#	    }
#	    if($j<($neventstot-1)){
#		$j++;
#		print KS ", $sorted_data[$j]";
#	    }
#	    print KS "\n";
#	}
#    }
#    print KS "     +/\n";
    print KS "      IF (V(1).le.EVENTS(1)) THEN\n";
    print KS "         ks = 0.\n";
    print KS "      ELSE\n";
    print KS "      DO 10 I = 1, NPAR\n";
    print KS "         IF (V(1).gt.EVENTS(I)) THEN\n";
#    print KS "            ks = (I-1)*1.\n";
    print KS "            ks = (I)*1.\n";
    print KS "         ELSE\n";
    print KS "            goto 15\n";
    print KS "         END IF\n";
    print KS "   10 CONTINUE\n";
    print KS "   15 END IF\n";
    print KS "      ks = ks/NPAR\n";
    print KS "      RETURN \n";
    print KS "      END\n";
    close KS;

#Compile and evaluate the shape and the cumulative distribution functions

#system "f77 -c $outfile";
#system "f77 -c shape.F";
`f77 -g -o shape $outfile shape.F`;
`shape > input.dat`;

`f77 -g -c ks.F`;
`f77 -g -c ksshape.F`;
`f77 -g -o ksshape ksshape.o ks.o`;
`ksshape > ksinput.dat`;



######################################################
#Step 4.5 - Make the final function
#		a) This is Linear interpolation with constant bins
#               b) The Y values of the function are writen in a data / / 
#                  statement in the function itself.

$evaluatedfunction = "input.dat";
open INTERPOLATE, "$evaluatedfunction";

while(<INTERPOLATE>) {
chomp($yi = $_);
push @y, $yi;
}
close INTERPOLATE;

$thefunction = "$func[$fn]";
open THEFUNCTION, ">$thefunction";

print THEFUNCTION "      DOUBLE PRECISION FUNCTION $funcname[$fn](V)\n";
print THEFUNCTION "      DOUBLE PRECISION V\n";
print THEFUNCTION "      DOUBLE PRECISION Yx, LBOUND, UBOUND, RANGE\n";
$evaluatedpoints = scalar(@y);
print THEFUNCTION "      DOUBLE PRECISION Y($evaluatedpoints )\n";
print THEFUNCTION "      INTEGER NBINS, I\n";
$numberofbins = $evaluatedpoints - 1;
print THEFUNCTION "      DATA NBINS  / $numberofbins /\n";
print THEFUNCTION "      DATA LBOUND, UBOUND / $lbnd[$fn], $ubnd[$fn]/\n";
print THEFUNCTION "      DATA y /";
$j=0;
for ($j=0; $j<$evaluatedpoints ; $j++){
	if($j==0){
		print THEFUNCTION "$y[$j]";
		if($j<($evaluatedpoints -1)){
		    $j++;
		    print THEFUNCTION ", $y[$j]";
		}
		print THEFUNCTION "\n";
	}
	else{
		print THEFUNCTION "     +, $y[$j]";
		if($j<($evaluatedpoints -1)){
		    $j++;
		    print THEFUNCTION ", $y[$j]";
		}
		if($j<($evaluatedpoints -1)){
		    $j++;
		    print THEFUNCTION ", $y[$j]";
		}
		print THEFUNCTION "\n";
	}
} 
print THEFUNCTION "     +/\n";
print THEFUNCTION "      If((V.lt.Lbound).or.(V.ge.UBound)) THEN\n";
print THEFUNCTION "      Yx = 0.0\n";
print THEFUNCTION "      ELSE\n";
print THEFUNCTION "      RANGE = UBOUND - LBOUND\n";
print THEFUNCTION "      i = 1 + INT(NBINS*(V - LBOUND)/RANGE)\n";
print THEFUNCTION "      Xi = (i-1)*RANGE/NBINS + LBOUND\n";
print THEFUNCTION "      Yx = (y(i+1) - y(i))*NBINS*(V - Xi)/RANGE + y(i)\n";
print THEFUNCTION "      END IF\n";
    use Env qw(QNORM);
    if ($QNORM == 1) {
	print THEFUNCTION "      $funcname[$fn] = Yx/dble($binwidth)\n";
    }
    else {
	print THEFUNCTION "      $funcname[$fn] = Yx/$normforhist\n";
    }
print THEFUNCTION "      RETURN\n";
print THEFUNCTION "      END\n";



######################################################
#Step 5 - Make kumac to:
#		a) turn keys function to multiquadric
#		b) write PostScript File of Shapes
#		c) output multiquadric function
#		d) run rename.script for input to CLFFT
#		e) calculate Kolmogorov-Smirnov Test (Dobs & P(D>Dobs))
print LOGFILE "Norm is $normforhist\n";
$final = "final.kumac";

open FINAL, ">$final";
print FINAL "1dh 2000 'Shape to fit' $bins[$fn] $lbnd[$fn] $ubnd[$fn]\n"; 
foreach $thiscut(@cuts) {
    print FINAL $thiscut."\n";
}

for ($fi = 1; $fi < $nfiles[$fn]+1; $fi ++) {
    $fi2 = $fi -1;
    print FINAL "hi/file 1 $ntup[$fn][$fi2]\n";
    $hid = 1000 + $fi;
    print FINAL "1dh $hid '$ntup[$fn][$fi2] shape' $bins[$fn] $lbnd[$fn] $ubnd[$fn]\n";
    print FINAL "nt/pl $ntid[$fn][$fi2]".".$vari[$fn][$fi2] $self[$fn][$fi2] -$hid\n";
    print FINAL "close 1\n";
    $factor = $weights[$fn][$fi2]*$lumi[$fn]*$NORMSUM[$fn];
    print FINAL "HI/OPE/ADD $hid 2000 2000 $factor 1\n";
}


print FINAL "opt nsta\n";

#print FINAL "fun1 10 ks.for 300 $lbnd[$fn] $ubnd[$fn]\n";
print FINAL "sigma s=array(301)\n";
print FINAL "vec/read s 'ksinput.dat'\n";
#print FINAL "hi/get_vec/cont 10 s\n";
#print FINAL "sigma x=array(301)\n";
#print FINAL "hi/get_vec/absc 10 x\n";
print FINAL "sigma x=array(301,$lbnd[$fn]"."#$ubnd[$fn]".")\n";
#print FINAL "fun1 20 f.for 301 $lbnd[$fn] $ubnd[$fn]\n";
print FINAL "sigma f=array(301)\n";
#print FINAL "hi/get_vec/cont 20 f\n";
print FINAL "vec/read f 'input.dat'\n";
$temp = ($ubnd[$fn] - $lbnd[$fn]);
print FINAL "sigma p=quad(f,$temp"."/300.0)\n";
print FINAL "sigma d = array(301)\n";
print FINAL "sigma d = p/p(301)-s\n";
print FINAL "sigma pp = p/p(301)\n";
print FINAL "vec/write d 'difference.txt'\n";
print FINAL "1dh 30 'P(x) and S(x)' 300 $lbnd[$fn] $ubnd[$fn]\n";
print FINAL "1dh 31 'S(x)' 300 $lbnd[$fn] $ubnd[$fn]\n";
print FINAL "1dh 32 'f(x)' 300 $lbnd[$fn] $ubnd[$fn]\n";
print FINAL "1dh 33 'P(x)-S(x)' 300 $lbnd[$fn] $ubnd[$fn]\n";
print FINAL "hi/put_vec/cont 30 pp\n";
print FINAL "hi/put_vec/cont 31 s\n";
print FINAL "hi/put_vec/cont 32 f\n";
print FINAL "hi/put_vec/cont 33 d\n";
print FINAL "zone 2 2\n";
print FINAL "close 66\n";
($psfile,@bunk3) = split /\./, $func[$fn];
print FINAL "FORTRAN/FILE 66 $funcname[$fn]".".ps\n";
print FINAL "META 66 -111\n";
@name = split /\_/, $funcname[$fn];
$title = "";
foreach $part (@name) {
    $title = $title."$part ";
}
print FINAL "title '$title '\n";
print FINAL "hi/pl 30\n";
print FINAL "hi/pl 31 S\n";
print FINAL "hi/pl 33\n";
print FINAL "zone 1 2 2 s\n";
print FINAL "opt sta\n";
print FINAL "set dmod 1\n";
print FINAL "hi/pl 2000\n";
print FINAL "hi/pl 32 CS\n";
print FINAL "close 66\n";
#print FINAL "exit\n";
close FINAL;


`paw -b final > delme.txt`;
#system "rename.script $func[$fn]";
#system "cp -f f.for $func[$fn]";

$diftemp = `cat difference.txt`;
#$f0 = "f$fn".".for";
#system "cp f.for $f0";
unlink ("difference.txt", "ks.for") or warn "Can't clean up \n";
@diftemp2 = split /\s+/, $diftemp;
foreach $event(@diftemp2){
    if ($event =~ /E/) {
	($number,$exponent) = split /E/, $event;
	push @diftemp3, abs($number)*10**$exponent;
#		print abs($number)*10**$exponent."\n";
    }	
    else {
	push @diftemp3, abs($event);
#	    print abs($event)."\n";
    }
}

##########################################
#  Clean up

undef @diftemp;
undef @diftemp2;
@diftemp4 = sort numerically @diftemp3;
$dif = $diftemp4[scalar(@diftemp4)-1];
undef @diftemp3;
undef @diftemp4;
print LOGFILE "Dobs = $dif\n";
if ($neventstot != 0) {$rootNeff = sqrt($neventstot);}
else {$rootNeff = 1;}
$Prob = &Qks(($rootNeff + .12 + .11/$rootNeff)*$dif);
print LOGFILE "P(D > Dobs) = $Prob\n";
$temp = ($rootNeff + .12 + .11/$rootNeff)*$dif;
system "cp shape.F shape$func[$fn]";
unlink ("f.for", "final.kumac", "$funcname[$fn].dat") or warn "Can't clean up\n";
unlink ("shape.F", "shape", "$funcname[$fn].o", "shape.o", "input.dat");
unlink ("ksshape.F", "ksshape", "ksshape.o", "ks.F", "ks.o", "ksinput.dat");
unlink ("ks$funcname[$fn].dat", "$outfile", "keys$funcname.o");
unlink ("shape$func[$fn]", "keys$func[$fn]");

undef @data;
undef @temp;
undef @y;

#########################################
# Print a summary txt file

($sum,$bunk) = split /\./, $func[$fn];
$summ = "$sum".".txt";
open SUMMARY, ">$summ";
print SUMMARY "$func[$fn]\n";
print SUMMARY "Total number of events            = $neventstot\n";
for($i = 1; $i<$nfiles; $i++) {
    print SUMMARY "Events in File $i                 = $nevents[$fn]\n";
}
print SUMMARY "Dobs                              = $dif\n";
print SUMMARY "Lamda                             = $temp\n";
print SUMMARY "Kolmogorov-Smirnov Test P(D>Dobs) = $Prob\n";
print SUMMARY "The function is normalized to fit a histogram with:\n";
print SUMMARY "          $bins[$fn] bins\n";
print SUMMARY "          $lbnd[$fn] lower bound\n";
print SUMMARY "          $ubnd[$fn] upper bound\n";
print SUMMARY "          $NPAReff Effective number of events (used for KS prob)\n";
$eventsexpected = $normforhist/($ubnd[$fn]-$lbnd[$fn])*$bins[$fn];
print SUMMARY "The expected number of events is  $eventsexpected\n";
#for ($ii = 0; $ii < $nfiles[$fn]; $ii++) {
#}

close SUMMARY;
}
#	##############################################
#		Task 3 - put all into 1 big file for CLFFT
#	##############################################


#	##############################################
#		Task 2 - Make CLFFT user function
#	##############################################

close LOGFILE;
exit;

#this function is to calculate the P(D>Dobs) for the KS test
sub Qks {
    $lambda = shift;
    $Q=0;
    for ($i = 1; $i<1000; $i++) {
	$Q += ((-1.)**($i-1.))*exp(-2.*($i*$lambda)**2.);
    }
    $Q *= 2.;
    return $Q;
}
