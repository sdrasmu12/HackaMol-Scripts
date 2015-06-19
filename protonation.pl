use HackaMol;
use HackaMol::X::NERF;
use Modern::Perl;

my $mol = HackaMol->read_file_mol('cystine_hvy.xyz');
$mol->print_xyz;1

#build in two oxygen
my @atoms = map{$_->xyz} $mol->all_atoms;
my $nerf = HackaMol::X::NERF->new;

my @atoms = map{$_->xyz} $mol->all_atoms;
my $nerf = HackaMol::X::NERF->new;

my $o1 = $nerf->extend_abc( $atoms[3], $atoms[1], $atoms[2], 1.3, 120, 180 );
my $o2 = $nerf->extend_abc( $atoms[9], $atoms[7], $atoms[8], 1.3, 120, 180 );

$mol->push_atoms (
    map{ HackaMol::Atom->new(Z => 8, coords => [$_]) } ($o1,$o2)
);
$mol->print_xyz;

my @atoms = map{$_->xyz} $mol->all_atoms;
my $nerf = HackaMol::X::NERF->new;

my $h1 = $nerf->extend_abc( $atoms[2], $atoms[1], $atoms[0], 1.0, 120, 1 );
my $h2 = $nerf->extend_abc( $atoms[2], $atoms[1], $atoms[0], 1.0, 120, 180 );
my $h3 = $nerf->extend_abc( $atoms[7], $atoms[5], $atoms[6], 1.0, 120, 1 );
my $h4 = $nerf->extend_abc( $atoms[7], $atoms[5], $atoms[6], 1.0, 120, 180 );
my $h5 = $nerf->extend_abc( $atoms[2], $atoms[0], $atoms[1], 1.47, 109, 120);
my $h6 = $nerf->extend_abc( $atoms[8], $atoms[6], $atoms[7], 1.47, 109, 90);
my $h7 = $nerf->extend_abc( $atoms[5], $atoms[3], $atoms[4], 1.47, 109, -90);
my $h8 = $nerf->extend_abc( $atoms[5], $atoms[3], $atoms[4], 1.47, 109, 90);
my $h9 = $nerf->extend_abc( $atoms[11], $atoms[9], $atoms[10], 1.47, 109, 270);
my $h10 = $nerf->extend_abc( $atoms[11], $atoms[9], $atoms[10], 1.47, 109, 90);
my $h11 = $nerf->extend_abc( $atoms[3], $atoms[2], $atoms[12], 0.96, 180, 180);
my $h12 = $nerf->extend_abc( $atoms[9], $atoms[8], $atoms[13], 0.96, 180, 90);


$mol->push_atoms (
    map{ HackaMol::Atom->new(Z => 1, coords => [$_]) } ($h1,$h2,$h3,$h4,$h5,$h6,$h7,$h8,$h9,$h10,$h11,$h12)
);
$mol->print_xyz;

$mol->print_xyz('cys_oo.xyz');1;

