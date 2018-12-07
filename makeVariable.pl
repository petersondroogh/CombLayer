#!/usr/bin/perl -w

use lib "/home/ansell/exe/C++";
use strict;
use CVar;

sub filterPrivate
{
  my $Fname=shift;
  my @nameVec;
  my @matVec;

  open(DX,$Fname) or die ("No file:".$Fname);

  my $line;
  my $pVal=0;
  while(defined($line=<DX>) && $pVal<2) 
  {
    $pVal=1 if ($line=~/private:/);
    $pVal=2 if ($line=~/public:/);
    if ($line=~/double\s+([a-zA-Z]+)\;/)
    {
      push(@nameVec,$1);
    }
    if ($line=~/int\s+([a-zA-Z]+)\;/)
    {
      my $M=$1;
      if ($M=~/[mM]at/)
        {
	  push(@matVec,$M);
	}
    }
  }
  return (\@nameVec,\@matVec);
}



my ($nameVec,$matVec)=filterPrivate($ARGV[0]);


foreach my $item (@$nameVec)
{
  print lcfirst($item),"=Control.EvalVar<double>(keyName+\"",
      ucfirst($item),"\");\n";
}

foreach my $item (@$matVec)
{
  print lcfirst($item),
      "=ModelSupport::EvalMat<int>(Control,keyName+\"",
      ucfirst($item),"\");\n";
}

