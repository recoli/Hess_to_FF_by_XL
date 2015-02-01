#!/usr/bin/perl

#*************************************************************************
#                ==== HESSIAN >>> to >> FF ====
#  Convert Hessian to bond stretching and angle bending parameters
#  Hessian is read from Gaussian fchk file
#  Output parameters are in GROMACS format
#
#  Copyright (C) 2012-2015 Xin Li <lixin.reco@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#*************************************************************************

#  Usage: perl Hess_to_FF_by_XL.pl xxx.fchk [2-3 atoms]
#  Ref: JM Seminario, Int. J. Quantum Chem., 1996, 60, 1271-1277.

use strict;
use warnings;
use Math::Trig;

if (@ARGV==0) { exit; }

# Gaussian fchk file

my $fchk = shift(@ARGV);
die "File $fchk doesn't exist!\n" unless (-e $fchk);

# get atom index and job type (bond, angle, dihedral)

die "Please specify 2-3 atoms!\n" if (@ARGV<2 or @ARGV>3);

my ($ai,$aj,$ak,$al,$jobtype);

$ai = $ARGV[0] - 1;
$aj = $ARGV[1] - 1;
$jobtype = "BOND";

if (@ARGV==3)
{
   $ak = $ARGV[2] - 1;
   $jobtype = "ANGLE";
}

# read fchk file for coordinates and Hessian

my (@g09hess, @hess, @g09coord, @coord);
my ($i, $j, $natoms, $ndf);

open FL,"$fchk" or die;
while(<FL>)
{
   # read coordinates
   if (/Current cartesian coordinates/)
   {
      # number of atoms
      $ndf = (split)[5];
      $natoms = $ndf/3;

      # 3N coordinates
      $i = 0;
      while(<FL>)
      {
         @_ = (split);
         for $j ( 0 .. @_-1 )
         {
            $g09coord[$i] = $_[$j];
            $i++;
         }
         last if $i==$ndf;
      }
      die "Error in reading coordiantes!\n" unless $i==$ndf;
   }

   # read hessian, 3N*(3N+1)/2 elements
   if (/Cartesian Force Constants/)
   {
      $i = 0;
      while (<FL>)
      {
         @_ = (split);
         for $j ( 0 .. @_-1 )
         {
            $g09hess[$i] = $_[$j];
            $i++;
         }
         last if $i==$ndf*($ndf+1)/2;
      }
      die "Error in reading Hessian!\n" unless $i==$ndf*($ndf+1)/2;
   }
}
close FL;

# get coordinates

for $i ( 0 .. $natoms-1 )
{
   $coord[$i][0] = $g09coord[$i*3];
   $coord[$i][1] = $g09coord[$i*3+1];
   $coord[$i][2] = $g09coord[$i*3+2];
}

# compute vectors

my (@vec_ij, @vec_ji, $r_ij, $r_ji);
my (@vec_jk, @vec_kj, $r_jk, $r_kj);
my (@vec_kl, @vec_lk, $r_kl, $r_lk);
my (@vec_PA, @vec_PC, @vec_NABC, @vec_NBCD);

# i->j
$vec_ij[0] = $coord[$aj][0] - $coord[$ai][0];
$vec_ij[1] = $coord[$aj][1] - $coord[$ai][1];
$vec_ij[2] = $coord[$aj][2] - $coord[$ai][2];
$r_ij = dist(@vec_ij);
@vec_ij = norm(@vec_ij);

my ($r0, $a0);
$r0 = $r_ij;
$r0 *= (0.52917721092 * 0.1); # nm

if ($jobtype eq "ANGLE")
{
   # k->j
   $vec_kj[0] = $coord[$aj][0] - $coord[$ak][0];
   $vec_kj[1] = $coord[$aj][1] - $coord[$ak][1];
   $vec_kj[2] = $coord[$aj][2] - $coord[$ak][2];
   $r_kj = dist(@vec_kj);
   @vec_kj = norm(@vec_kj);

   $a0 = acos(inn_p(@vec_ij, @vec_kj))/pi * 180.0;
   
   my @vec_N = norm(out_p(@vec_kj, @vec_ij));
   @vec_PA = out_p(@vec_N, @vec_ij);
   @vec_PC = out_p(@vec_kj, @vec_N);
}

# get full Hessian

my $count = 0;
for $i ( 0 .. $ndf )
{
   for $j ( 0 .. $i )
   {
      $hess[$i][$j] = $g09hess[$count];
      $hess[$j][$i] = $hess[$i][$j] if $i!=$j;
      $count++;
   }
}

# Calculate hessian matrix for pair potential

my (@AB, @CB);

for $i ( 0 .. 2 )
{
   for $j ( 0 .. 2 )
   {
      $AB[$i*3+$j] = -$hess[$ai*3+$i][$aj*3+$j];
      $CB[$i*3+$j] = -$hess[$ak*3+$i][$aj*3+$j] if ($jobtype eq "ANGLE");
   }
}

# calculate force constant

if ($jobtype eq "BOND")
{
   my $kbond = CalcBond(@AB, @vec_ij);
   printf "%-8d%-8d%-8d%15.6e%15.6e\n", 
          $ai+1, $aj+1, 1, $r0, $kbond;
}
elsif ($jobtype eq "ANGLE")
{
   my $kangle = CalcAngle(@AB, @CB, @vec_PA, @vec_PC);
   printf "%-7d%-7d%-7d%-7d%15.6e%15.6e\n", 
          $ai+1, $aj+1, $ak+1, 1, $a0, $kangle;
}

### subroutines

# calculate eigenvalue of a 3x3 matrix

sub GetEigValue
{
   unless (@_-9==0)
   {
      die "Error: Only 3x3 matrix is supported in GetEigValue!\n";
   }

   # save 9 elements in @A

   my ($i,@A);
   for $i ( 0 .. 8 )
   {
      $A[$i] = $_[$i];
   }

   # compute the eigenvalues
   # http://en.wikipedia.org/wiki/Eigenvalue_algorithm
    
   my $q = m_tr(@A)/3.0;
   my @AqI;
   for $i ( 0 .. 8 )
   {
      $AqI[$i] = $A[$i];
      $AqI[$i] -= $q if ($i==0 or $i==4 or $i==8);
   }
   my $p = m_tr(m_mult(@AqI,@AqI));
   $p = sqrt($p/6.0);

   my @B;
   for $i ( 0 .. 8 )
   {
      my $tmp = 0.0;
      $tmp = 1.0 if ($i==0 or $i==4 or $i==8);
      $B[$i] = 1.0/$p * ($A[$i]-$q*$tmp);
   }

   my $r = calcDET(@B)/2.0;
   if ($r<=-1)
   {
      #printf "Warning: \$r=$r is smaller than -1.0!\n";
      $r = -1;
   }
   elsif ($r>=1)
   {
      #printf "Warning: \$r=$r is greater than 1.0!\n";
      $r = 1;
   }
   my $phi = acos($r)/3;

   my ($eig0,$eig1,$eig2);
   $eig0 = $q + 2.0 * $p * cos($phi);
   $eig2 = $q + 2.0 * $p * cos($phi + pi * (2.0/3.0));
   $eig1 = 3.0 * $q - $eig0 - $eig2;

   # compute eigenvectors and print results

   my @eigvec = (GetEigVector(@A,$eig0),GetEigVector(@A,$eig1),GetEigVector(@A,$eig2));

   return ($eig0,$eig1,$eig2,@eigvec);
}

# calculate eigenvectors of a 3x3 matrix

sub GetEigVector
{
   # save 9 elements of (A-lambda*I) in @tmp
   
   my ($i,@tmp);
   for $i ( 0 .. 8 )
   {
      $tmp[$i] = $_[$i];
      $tmp[$i] -= $_[9] if ($i==0 or $i==4 or $i==8);
   }
   
   # calculate eigenvector
   
   my (@a,@b1,@b2);
   
   $a[0] = $tmp[0];
   $a[1] = $tmp[1];
   $a[2] = $tmp[3];
   $a[3] = $tmp[4];
   
   $b1[0] = -$tmp[2];
   $b1[1] =  $tmp[1];
   $b1[2] = -$tmp[5];
   $b1[3] =  $tmp[4];
   
   $b2[0] =  $tmp[0];
   $b2[1] = -$tmp[2];
   $b2[2] =  $tmp[3];
   $b2[3] = -$tmp[5];
   
   my $x = calcDET(@b1)/calcDET(@a);
   my $y = calcDET(@b2)/calcDET(@a);
   my $z = 1.0;
   
   # check the result
   my $check = $tmp[6]*$x + $tmp[7]*$y + $tmp[8]*$z;
   if (abs($check)>1.0e-3)
   {
      printf "Error: \$check in subroutine getEigVec is not zero! %e\n", $check;
      exit;
   }
   
   my $r = sqrt($x**2+$y**2+$z**2);
   $x /= $r;
   $y /= $r;
   $z /= $r;
   
   return ($x,$y,$z);
}

# calculate the determinant of a 2x2 or 3x3 matrix

sub calcDET
{
   my $r;
   if (@_-4==0)
   {
      $r = $_[0]*$_[3] - $_[1]*$_[2];
   }
   elsif (@_-9==0)
   {
      $r = $_[0]*$_[4]*$_[8] + $_[1]*$_[5]*$_[6] + $_[2]*$_[3]*$_[7]
         - $_[2]*$_[4]*$_[6] - $_[1]*$_[3]*$_[8] - $_[0]*$_[5]*$_[7];
   }
   else
   {
      die "Error: Only 2x2 and 3x3 matrices are supported in subroutine calcDET!\n";
   }
   return $r;
}

# calculate the trace of a 3x3 matrix

sub m_tr
{
   die "Incorrect number of elements in subroutine m_tr!\n" 
      unless @_==9;
   my $c = $_[0]+$_[4]+$_[8];
   return ($c);
}

# multiplication of 3x3 matrices

sub m_mult
{
   die "Incorrect number of elements in subroutine mmult!\n"
      unless @_==18;
   my ($i,$j,$k);
   my (@a,@b,@c);
   for $i ( 0 .. 8 )
   {
      $a[$i] = $_[$i];
      $b[$i] = $_[$i+9];
   }
   for $i ( 0 .. 2 )
   {
      for $j ( 0 .. 2 )
      {
         $c[$i*3+$j] = 0.0;
         for $k ( 0 .. 2 )
         {
            $c[$i*3+$j] += $a[$i*3+$k] * $b[$k*3+$j];
         }
      }
   }
   return (@c);
}

# inner product of two vectors

sub inn_p
{
   my ($i);
   my (@a,@b,$c);
   for $i ( 0 .. 2 )
   {
      $a[$i] = $_[$i];
      $b[$i] = $_[$i+3];
   }
   $c = $a[0]*$b[0] + $a[1]*$b[1] + $a[2]*$b[2];
   return ($c);
}

# outer product of two vectors

sub out_p
{
   my ($i);
   my (@a,@b,@c);
   for $i ( 0 .. 2 )
   {
      $a[$i] = $_[$i];
      $b[$i] = $_[$i+3];
   }
   $c[0] = $a[1]*$b[2] - $a[2]*$b[1];
   $c[1] = $a[2]*$b[0] - $a[0]*$b[2];
   $c[2] = $a[0]*$b[1] - $a[1]*$b[0];
   return (@c);
}

# normalization of a vector

sub norm
{
   my @a;
   my $r = sqrt($_[0]**2+$_[1]**2+$_[2]**2);
   for $i ( 0 .. 2 )
   {
      $a[$i] = $_[$i]/$r;
   }
   return (@a);
}

# length of a vector

sub dist
{
   my $r = sqrt($_[0]**2+$_[1]**2+$_[2]**2);
   return ($r);
}

# calculate bond stretching parameters

sub CalcBond
{
   my ($i,@AB,@vec_ij);
   for $i ( 0 .. 8 )   {   $AB[$i] = shift(@_);   }
   for $i ( 0 .. 2 )   {   $vec_ij[$i] = shift(@_);   }

   # compute the eigenvalues of @AB
    
   @_ = GetEigValue(@AB);
   my $eig0 = shift(@_);
   my $eig1 = shift(@_);
   my $eig2 = shift(@_);
   my @eigvec = (@_);
   
   my $k = $eig0 * abs($vec_ij[0]*$eigvec[0] + $vec_ij[1]*$eigvec[1] + $vec_ij[2]*$eigvec[2])
         + $eig1 * abs($vec_ij[0]*$eigvec[3] + $vec_ij[1]*$eigvec[4] + $vec_ij[2]*$eigvec[5])
         + $eig2 * abs($vec_ij[0]*$eigvec[6] + $vec_ij[1]*$eigvec[7] + $vec_ij[2]*$eigvec[8]);
   
   $k *= 6.02214129 * 4.35974434 * 100.0; # kJ/mol Bohr^-2
   $k /= (0.52917721092 * 0.1) ** 2;      # kJ/mol nm^-2
   return ($k);
}

# calculate angle bending parameters

sub CalcAngle
{
   my ($i, @AB, @CB, @vec_PA, @vec_PC);
   for $i (0 .. 8) { $AB[$i] = shift(@_); }
   for $i (0 .. 8) { $CB[$i] = shift(@_); }
   for $i (0 .. 2) { $vec_PA[$i] = shift(@_); }
   for $i (0 .. 2) { $vec_PC[$i] = shift(@_); }

   @_ = GetEigValue(@AB);
   my $eig0 = shift(@_);
   my $eig1 = shift(@_);
   my $eig2 = shift(@_);
   my @eigvec = (@_);
   my $k1 = $eig0 * abs($vec_PA[0]*$eigvec[0] + $vec_PA[1]*$eigvec[1] + $vec_PA[2]*$eigvec[2])
          + $eig1 * abs($vec_PA[0]*$eigvec[3] + $vec_PA[1]*$eigvec[4] + $vec_PA[2]*$eigvec[5])
          + $eig2 * abs($vec_PA[0]*$eigvec[6] + $vec_PA[1]*$eigvec[7] + $vec_PA[2]*$eigvec[8]);

   @_ = GetEigValue(@CB);
   $eig0 = shift(@_);
   $eig1 = shift(@_);
   $eig2 = shift(@_);
   @eigvec = (@_);
   my $k2 = $eig0 * abs($vec_PC[0]*$eigvec[0] + $vec_PC[1]*$eigvec[1] + $vec_PC[2]*$eigvec[2])
          + $eig1 * abs($vec_PC[0]*$eigvec[3] + $vec_PC[1]*$eigvec[4] + $vec_PC[2]*$eigvec[5])
          + $eig2 * abs($vec_PC[0]*$eigvec[6] + $vec_PC[1]*$eigvec[7] + $vec_PC[2]*$eigvec[8]);
   
   my $kangle = 1.0 / (1.0 / ($r_ij**2 * $k1) + 1.0 / ($r_kj**2 * $k2));
   $kangle *= 6.02214129 * 4.35974434 * 100.0;      # kJ/mol
   return ($kangle);
}

