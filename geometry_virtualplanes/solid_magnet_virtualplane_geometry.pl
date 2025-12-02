#!/usr/bin/perl -w
# use strict;
use warnings;
our %detector;
our %configuration;
our %parameters;

use Getopt::Long;
use Math::Trig;

my $DetectorName = 'solid_magnet_virtualplane';

# my $DetectorMother="root";
my $DetectorMother="cc_pro_tcd";

my $x1	= 0;
my $x2  = 0;
my $x3  = 0;
my $x4  = -10;
my $x5  = 10;
my $y4  = 0;
my $y5  = 0;
my $z1	= -98;
my $z3  = -167;
my $z2  = -131;
my $z4  = -114.5;
my $z5  = -114.5;
my $z6  = -267;
my $hx	= 50;
my $hy	= 50;
my $hx4	= 16.4;
my $hy4	= 16.4;
my $hx5	= 16.4;
my $hy5	= 16.4;

sub solid_magnet_virtualplane
{
make_1();
make_2();
make_3();
make_4();
make_5();
make_6();
}

sub make_1
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_1";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x1*cm 0*cm $z1*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx*cm $hy*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 51;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}

sub make_2
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_2";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x2*cm 0*cm $z2*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx*cm $hy*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 52;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_3
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_3";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x3*cm 0*cm $z3*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx*cm $hy*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 53;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_4
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_4";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x4*cm $y4*cm $z4*cm";
 $detector{"rotation"}    = "0*deg 90*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx4*cm $hy4*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 54;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_5
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_5";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x5*cm $y5*cm $z5*cm";
 $detector{"rotation"}    = "0*deg 90*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx5*cm $hy5*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 55;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}
sub make_6
{
 my %detector=init_det();
 $detector{"name"}        = "$DetectorName\_6";
 $detector{"mother"}      = "$DetectorMother";
 $detector{"description"} = $detector{"name"};
 $detector{"pos"}         = "$x3*cm 0*cm $z6*cm";
 $detector{"rotation"}    = "0*deg 0*deg 0*deg";
 $detector{"color"}       = "CC6633";
 $detector{"type"}       = "Box";
 $detector{"dimensions"} = "$hx*cm $hy*cm 0.0001*cm";	    
 $detector{"material"}    = "G4_Galactic";
 $detector{"mfield"}      = "no";
 $detector{"ncopy"}       = 1;
 $detector{"pMany"}       = 1;
 $detector{"exist"}       = 1;
 $detector{"visible"}     = 1;
 $detector{"style"}       = 0;
 $detector{"sensitivity"} = "flux";
 $detector{"hit_type"}    = "flux";
 my $ID = 56;
 $detector{"identifiers"} = "id manual $ID";
 print_det(\%configuration, \%detector);
}


solid_magnet_virtualplane();
1;
