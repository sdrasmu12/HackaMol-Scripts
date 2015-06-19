use Modern::Perl;
use HackaMol::X::Orca;
use Math::Vector::Real;
use HackaMol;
use Time::HiRes qw(time);
use Data::Dumper;

my $t1 = time;


my $mol = HackaMol->new->read_file_mol(shift);

$mol->push_charges(0);
$mol->multiplicity(1);

my $dihe = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(3),
                                $mol->get_atoms(2),
                                $mol->get_atoms(1), 
                                $mol->get_atoms(0),
                                             ],
);


#my $bond = HackaMol::Bond->new(
#                                is_constrained => 1,
#                                    atoms=> [
#                                $mol->get_atoms(2),
#                                $mol->get_atoms(1),
#                                             ],
#);

#$mol->push_bonds($bond);
$mol->push_dihedrals($dihe);
#my @dihes = HackaMol->new->build_dihedrals ($mol->select_atoms( sub{ $_->Z != 1 } ) );
#$_->is_constrained(1) foreach @dihes;
#$mol->push_dihedrals(@dihes);


my $orca2 = HackaMol::X::Orca->new(
      mol             => $mol,
#      has_constraints => 1,      
      theory          => 'HF-3c',
      exe             => '/Users/chem_student/perl5/apps/orca_3_0_3_macosx_openmpi165/orca',
      scratch         => 'tmp',
);

$mol->gt(0);
my $fh = $mol->print_xyz('shit.xyz');
open(my $in, ">", "plotly.txt") or die "couldn't open";
my $min = 9999;
#unrestrained ones
foreach (0,1){
  $mol->gt($_);
  my @energies = $orca2->opt;
  #my @energies = $orca2->ener;
  $min = $energies[-1] if ($min > $energies[-1] );
  my $e_rel = ($energies[-1]-$min)*627.51;  
  printf ("%10.3f %10.2f\n", $dihe->dihe_deg, $energies[-1]*627.51);
  printf $in ("%10.3f %10.2f\n", $dihe->dihe_deg, $e_rel, $energies[-1]*627.51);  
  #my $mol2 = HackaMol->new->read_file_mol('tmp/mol.xyz');
  #$mol2->print_xyz;
  #$mol2->print_xyz($fh);
}

my $mol1 = HackaMol->new->read_file_mol('tmp/mol.xyz');
$mol1->print_xyz($fh);

#      has_constraints => 1,      

$orca2->has_constraints(1);

foreach (2 .. $mol->tmax){
  $mol->gt($_);
  my @energies = $orca2->opt;
  my $e_rel = ($energies[-1]-$min)*627.51;
  printf ("%10.3f %10.2f %14.6f\n", $dihe->dihe_deg, $e_rel, $energies[-1]*627.51);
  printf $in ("%10.3f %10.2f %14.6f\n", $dihe->dihe_deg, $e_rel, $energies[-1]*627.51);
  my $mol2 = HackaMol->new->read_file_mol("tmp/mol.xyz");
  $mol2->print_xyz($fh);
}

#$orca->map_input;
#$orca2->load_engrad;

my $t2 = time;

printf ("%10.2f\n", $t2-$t1);
