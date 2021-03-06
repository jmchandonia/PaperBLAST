#!/usr/bin/perl -w
# The last phase of the pipeline: building the database
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use pbutils;
use DBI;

my @allsteps = qw{sprot db stats compare};
my $allsteps = join(",",@allsteps);
my $comparedir = "$Bin/../data";

my $usage = <<END
run_final.pl -in downloads [ -work work ] [ -test ] [ -steps ... ]
	[ -compare $comparedir ]

Given the downloads and work directories, and the output of
run_snippets.pl, builds a database in work/data, and compares it to
the existing database (from $comparedir)

END
;

my ($test, $indir);
my $workdir = "work";
my $dosteps = $allsteps;

sub maybe_run($) {
    my ($cmd) = @_;
    if (defined $test) {
        print STDERR "Would run\t$cmd\n";
    } else {
        print STDERR "Running $cmd\n";
        system($cmd) == 0
            || die "Error running $cmd: $!\n";
    }
}

die $usage
    unless GetOptions('in=s' => \$indir,
                      'work=s' => \$workdir,
                      'steps=s' => \$dosteps,
                      'test' => \$test,
                      'compare=s' => \$comparedir)
    && @ARGV == 0
    && defined $indir && defined $workdir;
die "No such directory: $indir\n" unless -d $indir;
die "No such directory: $workdir\n" unless -d $workdir;
die "No such directory: $comparedir\n" unless -d $comparedir;

my @dosteps = split /,/, $dosteps;
my %dosteps = map { $_ => 1 } @dosteps;
my %allsteps = map { $_ => 1 } @allsteps;
foreach my $step (keys %dosteps) {
    die "Unrecognized step: $step\n" unless exists $allsteps{$step};
}

print STDERR "Test mode\n" if defined $test;

my $outdir = "$workdir/data";
&mkdir_if_needed($outdir);

if (exists $dosteps{sprot}) {
  my $sprot = "$indir/uniprot_sprot.dat.gz";
  die "No such file: $sprot\n" unless -e $sprot;
  &maybe_run("(zcat $sprot | $Bin/sprotCharacterized.pl > $workdir/sprot.char.tab) >& $workdir/sprot.char.log");
}

if (exists $dosteps{db}) {
  die "No such directory: $indir/ecocyc/data\n" unless -d "$indir/ecocyc/data";
  foreach my $file (qw{snippets_comb hits.papers pmclinks.papers generif_tab.rif}) {
    die "No such file: $workdir/$file" unless -e "$workdir/$file";
  }
  &maybe_run("$Bin/buildLitDb.pl -dir $outdir -snippets $workdir/snippets_comb -rif $workdir/generif_tab.rif -sprot $workdir/sprot.char.tab -ecocyc $indir/ecocyc/data $workdir/hits $workdir/pmclinks $workdir/generif_tab");
}

if (exists $dosteps{stats}) {
  if (!defined $test) {
    my $sqldb = "$outdir/litsearch.db";
    my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
    my $pm1 = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT pmId FROM GenePaper WHERE pmId <> "" });
    my $pm2 = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT pmId FROM GeneRIF });
    my $pmc = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT pmcId FROM GenePaper
                                           WHERE pmId = "" AND pmcId <> "" });
    my %pm = map { $_ =>  1 } @$pm1;
    foreach my $pm (@$pm2) { $pm{$pm} = 1; }
    my $nPaper = scalar(keys %pm) + scalar(@$pmc);

    my $nSeq = `grep -c ">" $outdir/uniq.faa`;
    chomp $nSeq;
    my $date = `date -r $workdir/epmc '+%B %-d %Y'`;
    chomp $date;
    my %stats = ( "nSeq" => $nSeq, "nPaper" => $nPaper, "date" => $date );
    open(STATS, ">", "$outdir/stats") || die "Error writing to $outdir/stats";
    while (my ($key,$value) = each %stats) {
      print STATS "$key\t$value\n";
    }
    close(STATS) || die "Error writing to $outdir/stats\n";
    $dbh->disconnect();
  }
  &maybe_run("cat $outdir/stats");
}

if (exists $dosteps{compare}) {
  my $newdb = "$outdir/litsearch.db";
  my $olddb = "$comparedir/litsearch.db";
  if (! -e $olddb) {
    print STDERR "No such file: $olddb -- skipping the comparison step\n";
  } else {
    &maybe_run("$Bin/compare_dbs.pl -old $olddb -new $newdb -out $workdir/dbcmp");
  }
}
