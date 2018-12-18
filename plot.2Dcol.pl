#!/usr/bin/perl
use Dislin;

# Set title
$ctit1 = 'Response function';
$twopi=2;

# File name
$file = $ARGV[0];
$min1 = $ARGV[1];
$max1 = $ARGV[2];
$min2 = $ARGV[3];
$max2 = $ARGV[4];
$step1 = $ARGV[5];
$step2 = $ARGV[6];
$levels = $ARGV[7];
$width = $ARGV[8];
$hight = $ARGV[9];

$xshift = 0;
$yshift = 0;
$zmin=0;
$zmax=0;
$level1=5.0;
if ($levels==18){
  @clines=qw/ -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 /;
} elsif ($levels==10){
  @clines=qw/ -0.2 -0.16 -0.12 -0.08 -0.04 0.04 0.08 0.12 0.16 0.2 /;
} elsif ($levels==20){
  @clines=qw/ -0.50 -0.45 -0.40 -0.35 -0.30 -0.25 -0.20 -0.15 -0.10 -0.05 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 /;
} elsif ($levels==16){
  @clines=qw/ -0.64 -0.32 -0.16 -0.08 -0.04 -0.02 -0.01 -0.005 0.005 0.01 0.02 0.04 0.08 0.16 0.32 0.64 /;
}


$cstep=256/($levels+1);

$ox=1234;
#$ox=12000;
$i=-1;
$n=0,$m=0;

print("Open file $file.\n");
open(FILEA, "<$file");
while (<FILEA>) {
  ($x,$y,$d,$z) = split;

  if ($x > $min1){  
     if ($x < $max1){
	  if ($y > $min2){
	      if ($y < $max2){
	    
		  if ($x != $ox){
		      $i++;
		      $ox=$x;
		      $j=0;
		  }
		  $xray[$i]=$x+$xshift;
		  $yray[$j]=$y+$yshift;
		  $zmat[$i][$j]=$z;

		  if ($z > $zmax) {$zmax=$z};
		  if ($z < $zmin) {$zmin=$z};
		  if ($i>$n) {$n=$i};
		  if ($j>$m) {$m=$j};
		  $j++;
	      }
	  }
      }
  }
}

print ("Read array with dimension: $n $m\n");
$m++,$n++;
Dislin::metafl ('eps');
Dislin::setfil ($file.'.eps');

Dislin::page(1700,1700);
Dislin::scrmod ('revers');
Dislin::disini ();
# Height sets textsize!!!
Dislin::height(70);
Dislin::complx ();
Dislin::linwid (6); 
Dislin::shdmod ('poly', 'contur');

Dislin::axspos (500, 1370);
#Dislin::axslen(1200,1200);
#Dislin::axslen($width,$hight);

Dislin::texmod ('ON');
Dislin::name   ('\huge $\omega_1/2\pi$c (cm$^{-1}$)', 'X');
Dislin::name   ('\huge $\omega_3/2\pi$c (cm$^{-1}$)', 'Y');
Dislin::name   ('Intensity', 'Z');

Dislin::intax  ();
Dislin::autres ($n, $m);
Dislin::ax3len (1100, 1100, 1100);
Dislin::ax3len ($width,$hight, 1100);
Dislin::labels('EXP','Z');
Dislin::labdig(2,'Z');
Dislin::setvlt ('SPEC');

if ($zmin<0){
    if ($zmax<-$zmin){
	$zmax=-$zmin;
    }
    if ($zmax>-$zmin){
	$zmin=-$zmax;
    }
}

Dislin::graf  ($min1+$xshift, $max1+$xshift, $min1+$xshift, $step1, $min2+$yshift, $max2+$yshift, $min2+$yshift, $step2);

$min=$min1,$max=$max1;
if ($min2>$min) {$min=$min2};
if ($max2<$max) {$max=$max2};
Dislin::rline ($min+$xshift,$min+$yshift,$max+$xshift,$max+$yshift);

Dislin::height (60);

for ($i = 0; $i < $levels; $i++) {
#  print("Drawing level $i\n");
#    $zlev = $zmin + ($i+0.5) * ($zmax-$zmin)/($levels);
  $zlev=$zmax*$clines[$i];
  $c=abs($clines[$i]/$clines[0]);
  $a=1;#(1-$c)/2;
  $b=(1-$c)*3/4;
#  $a=1;
#  $b=(1-$c);
  print("$a $b $c\n");
  Dislin::setrgb ($b,$b,$a);
  Dislin::solid;
  if ($zlev>0){
#    Dislin::dashl;
    Dislin::setrgb ($a,$b,$b);
  }

  print("$i $zlev\n");
  Dislin::labels ('NONE', 'CONTUR');
#  if ($zlev>$$min){
#    if ($lev<$zmax){
      printf("Drawing level $i!\n");
      Dislin::contur (\@xray, $n, \@yray, $m, \@zmat, $zlev);
#    }
#  }
  $zd[i]=$zlev;
}

Dislin::color  ('FORE');
Dislin::height (70);
Dislin::title  ();
Dislin::disfin ();

