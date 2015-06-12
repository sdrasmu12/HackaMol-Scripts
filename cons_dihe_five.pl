use HackaMol;
use HackaMol::X::Orca;
use Modern::Perl;
use POSIX;

# setup

my $t1 = time;

my $mol = HackaMol->new->read_file_mol(shift);

$mol->push_charges(0);
$mol->multiplicity(1);


#establish dihedrals

my $dihe1 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(1),
                                $mol->get_atoms(0),
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                             ],
);

my $dihe2 = HackaMol::Dihedral->new(
                               is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(0),
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                $mol->get_atoms(10),
                                             ],
);

my $dihe3 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                $mol->get_atoms(10),
                                $mol->get_atoms(9),
                                             ],
);

my $dihe1p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                   atoms=> [
                                $mol->get_atoms(7),
                                $mol->get_atoms(6),
                                $mol->get_atoms(9),
                                $mol->get_atoms(10),
                                             ],
);


my $dihe2p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(7),
                                $mol->get_atoms(9),
                                $mol->get_atoms(10),
                               $mol->get_atoms(4),
                                             ],
);



# push dihedrals
$mol->push_dihedrals($dihe1);
$mol->push_dihedrals($dihe2);
$mol->push_dihedrals($dihe3);
$mol->push_dihedrals($dihe1p);
$mol->push_dihedrals($dihe2p);

#checks
say $dihe1->dihe_deg;
say $dihe2->dihe_deg;
say $dihe1p->dihe_deg;
say $dihe2p->dihe_deg;
say $dihe3->dihe_deg;
#set up orca

my $orca2 = HackaMol::X::Orca->new(
      mol             => $mol,
      has_constraints => 1,      
      theory          => 'HF-3c',
      exe             => '/Users/chem_student/perl5/apps/orca_3_0_3_macosx_openmpi165/orca',
      scratch         => 'tmp',
);

# get ready for the foreach loop
$mol->gt(0);
my $fh = $mol->print_xyz('shit.xyz');
open(my $in, ">", "plotly.txt") or die "couldn't open";
my %results;
my $partit = 10;
my @energies = (0);

foreach (0 .. $mol->tmax){
  $mol->gt($_);
  my $chi1 = ceil($dihe1->dihe_deg/$partit);
  my $chi1p = ceil($dihe1p->dihe_deg/$partit);
  my $chi2 = ceil($dihe2->dihe_deg/$partit);
  my $chi2p = ceil($dihe2p->dihe_deg/$partit);
  my $chi3 = ceil($dihe3->dihe_deg/$partit);
  # sets variables to the appropriate angle rounded up and in partitions set above

  unless( exists( $results{$chi1}{$chi2}{$chi3}{$chi1p}{$chi2p} ) ) {

  $results{$chi1}{$chi2}{$chi3}{$chi1p}{$chi2p}++;
  
  my @energies = $orca2->opt;
  printf ("%10.3f %10.3f %10.3f %10.3f %10.3f %14.6f\n",$dihe1->dihe_deg, $dihe1p->dihe_deg, $dihe2->dihe_deg, $dihe2p->dihe_deg, $dihe3->dihe_deg, $energies[-1]*627.51);
  printf $in ("%10.3f %10.3f %10.3f %10.3f %10.3f %14.6f\n",$dihe1->dihe_deg, $dihe1p->dihe_deg, $dihe2->dihe_deg, $dihe2p->dihe_deg, $dihe3->dihe_deg, $energies[-1
]*627.51);  
  my $mol2 = HackaMol->new->read_file_mol("tmp/mol.xyz");
  $mol2->print_xyz($fh);
  
#  else(
#  $results{$chi1}{$chi2}{$chi3}{$chi1p}{$chi2p}++;
#  );
  }

};

my $t2 = time;

printf ("%10.2f\n", ($t2-$t1)/60);
