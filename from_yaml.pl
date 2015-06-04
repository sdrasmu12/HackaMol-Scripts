use Modern::Perl;
use HackaMol;

my $mol = HackaMol->new->read_file_mol(shift);

foreach my $t (0 .. $mol->tmax){
  $mol->gt($t);
  $mol->print_xyz;
}


