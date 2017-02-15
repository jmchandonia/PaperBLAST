#!/usr/bin/perl -w
#######################################################
## litSearch.cgi
##
## Copyright (c) 2017 University of California
##
## Authors: Morgan Price
#######################################################
#
# Optional CGI garameters:
# query (protein sequence in FASTA or UniProt or plain format)
# vimss -- specify the vimss id to search for (mostly to make testing easier)

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw{gettimeofday};
use DBI;

sub fail($);
sub simstring($$$$$$$$$$$);
sub SubjectToGene($);

my $tmpDir = "../tmp";
my $blastall = "../bin/blast/blastall";
my $nCPU = 4;
my $base = "../data";
my $blastdb = "$base/uniq.faa";
my $sqldb = "$base/litsearch.db";
die "No such executable: $blastall" unless -x $blastall;
die "No such file: $blastdb" unless -e $blastdb;
die "No such file: $sqldb" unless -e $sqldb;

my $cgi=CGI->new;
my $query = $cgi->param('query') || "";

my $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $mo_dbh = DBI->connect('DBI:mysql:genomics:pub.microbesonline.org', "guest", "guest")
    || die $DBI::errstr;

my $documentation = <<END
<H3><A NAME="works">How It Works</A></H3>

<P>PaperBLAST relies on a database of protein sequences that are linked to scientific publications. These links come from automated text searches against <A HREF="http://europepmc.org/">EuropePMC</A> and from the manually curated information in <A HREF="http://www.uniprot.org/">Swiss-Prot</A> (the manually reviewed part of UniProt). Given this database and a query sequence, we use <A HREF="https://en.wikipedia.org/wiki/BLAST">protein-protein BLAST</A> to find similar sequences with E &lt; 0.001. As of January 2017, PaperBLAST contains over 200,000 proteins and links to over 50,000 papers.

<P>To link proteins in sequenced genomes to papers in EuropePMC, we query every locus tag or <A HREF="https://www.ncbi.nlm.nih.gov/refseq/">RefSeq</A> protein id that appears in the open access part of EuropePMC. We obtain the protein sequences and identifiers from <A HREF="http://www.microbesonline.org/">MicrobesOnline</A> or from RefSeq. We use queries of the form "locus_tag AND genus_name" to try to ensure that the paper is actually discussing the gene as opposed to something else whose identifier happens to match a locus tag. Note that EuropePMC indexes secret papers as well as open-access papers, so some of the links may be to papers that you cannot read (and that it would be illegal for our computers to read).</P>
<P>We also query EuropePMC for every locus tag that appears in the 300 most-referenced genomes. So, a gene may appear in the PaperBLAST results even though all of the papers that mention it are secret.</P>

<P>Finally, we index proteins from Swiss-Prot if their curators identified experimental evidence for the protein\'s function, as indicated by evidence code ECO:0000269. Most of these entries include links to articles in <A HREF="http://www.ncbi.nlm.nih.gov/pubmed/">PubMed</A><sup>&reg;</sup>. (We also use PubMed<sup>&reg;</sup> abstracts to help identify snippets of text that contain the locus tag from the secret papers.)

<P>The code for PaperBLAST is available <A HREF="https://github.com/morgannprice/PaperBLAST">here</A>.

<H3><A NAME="secret">Secrets</A></H3>

<P>PaperBLAST cannot provide snippets for many of the papers that are published in non-open-access journals. This limitation applies even if the paper is marked as "free" on the publisher\'s web site and is available in PubmedCentral or EuropePMC. Because these papers are not open access, PaperBLAST cannot download these papers en masse and scan them for snippets. If a journal that you publish in is marked as "secret," please consider publishing elsewhere.

<center><A HREF="http://morgannprice.org/">Morgan Price</A><BR>
<A HREF="http://genomics.lbl.gov/">Arkin group</A><BR>
Lawrence Berkeley National Laboratory
</center>
END
    ;

my $title = "PaperBLAST";
# utf-8 because that is the encoding used by EuropePMC
print
    header(-charset => 'utf-8'),
    start_html($title);

my $seq;
my $def = "";

if ($cgi->param('vimss')) {
    my $locusId = $cgi->param('vimss');
    die "Invalid vimss argument" unless $locusId =~ m/^\d+$/;
    my ($aaseq) = $mo_dbh->selectrow_array(
        qq{SELECT sequence FROM Locus JOIN AASeq USING (locusId,version)
            WHERE locusId = ? AND priority=1 },
        {}, $locusId);
    die "Unknown VIMSS id $locusId" unless defined $aaseq;
    $query = ">VIMSS$locusId\n$aaseq\n";
}
my $hasDef = 0;
if ($query =~ m/[A-Za-z]/) {
    $seq = "";
    my @lines = split /[\r\n]+/, $query;
    $def = shift @lines if @lines > 0 && $lines[0] =~ m/^>/;
    $def =~ s/^>//;
    foreach (@lines) {
        s/[ \t]//g;
        s/^[0-9]+//; # leading digit/whitespace occurs in UniProt format
        next if $_ eq "//";
        &fail("Error: more than one sequence was entered.") if m/^>/;
        &fail("Unrecognized characters in $_") unless m/^[a-zA-Z*]*$/;
        s/[*]/X/g;
        $seq .= uc($_);
    }
    $hasDef = 1 if $def ne "";
    $def = length($seq) . " a.a." if $def eq "";
}

if (!defined $seq) {
    my $exampleId = "3615187";
    my $refseqId = "WP_012018426.1";
    print
        h2($title),
        p(b("Find papers about a protein or its homologs")),
        start_form( -name => 'input', -method => 'POST', -action => 'litSearch.cgi'),
        p("Enter a sequence in FASTA or Uniprot format: ",
          br(),
          textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
        p(submit('Search'), reset()),
        end_form,
        p("Example: an 'alcohol dehydrogenase'",
          small("(" .
                a({ -href => "https://www.ncbi.nlm.nih.gov/protein/$refseqId",
                    -style => "color: black;" },
                  "$refseqId" )
                . ")"),
          "is actually the regulator",
          a({ -href => "litSearch.cgi?vimss=$exampleId",
              -title => "Show PaperBLAST hits" },
            i("ercA")) . "."),
        $documentation;
} else {
    my $procId = $$;
    my $timestamp = int (gettimeofday() * 1000);
    my $filename = $procId . $timestamp;
    my $seqFile = "$tmpDir/$filename.fasta";
    my $seqlen = length($seq);

    my $initial = substr($seq, 0, 10);
    $initial .= "..." if $seqlen > 10;
    $initial = "$seqlen a.a., $initial" if $hasDef;

    print
        h2("PaperBLAST Hits for $def ($initial)"),
        "\n";
    open(SEQ, ">", $seqFile) || die "Cannot write to $seqFile";
    print SEQ ">$def\n$seq\n";
    close(SEQ) || die "Error writing to $seqFile";
    my $hitsFile = "$tmpDir/$filename.hits";
    system($blastall, "-p", "blastp", "-d", $blastdb, "-i", $seqFile, "-o", $hitsFile,
           "-e", 0.001, "-m", 8, "-a", $nCPU, "-F", "m S") == 0 || die "Error running blastall: $!";
    my @hits = ();
    open(HITS, "<", $hitsFile) || die "Cannot read $hitsFile";
    while(<HITS>) {
        chomp;
        my @F = split /\t/, $_;
        push @hits, \@F;
    }
    close(HITS) || die "Error reading $hitsFile";
    unlink($seqFile);
    unlink($hitsFile);
    my $nHits = scalar(@hits);
    if ($nHits == 0) {
        print p("Sorry, no hits to proteins in the literature.");
    } else {
        print p("Found $nHits similar proteins in the literature:"), "\n";

        my %seen_subject = ();
        foreach my $row (@hits) {
            my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,$queryStart,$queryEnd,$subjectStart,$subjectEnd,$eVal,$bitscore) = @$row;
            next if exists $seen_subject{$subjectId};
            $seen_subject{$subjectId} = 1;

            my $dups = $dbh->selectcol_arrayref("SELECT duplicate_id FROM SeqToDuplicate WHERE sequence_id = ?",
                                                {}, $subjectId);
            my @subject_ids = ($subjectId);
            push @subject_ids, @$dups;

            my @genes = map { &SubjectToGene($_) } @subject_ids;
            @genes = sort { $a->{priority} <=> $b->{priority} } @genes;
            my @headers = ();
            my @content = ();
            my %paperSeen = (); # to avoid redundant papers -- pmId.pmcId.doi => term => 1
            my %paperSeenNoSnippet = (); # to avoid redundant papers -- pmId.pmcId.doi => 1

            # Produce top-level and lower-level output for each gene (@headers, @content)
            # Suppress duplicate papers if no additional terms show up
            # (But, a paper could show up twice with two different terms, instead of the snippets
            # being merged...)
            foreach my $gene (@genes) {
                die "No subjectId" unless $gene->{subjectId};
                foreach my $field (qw{showName URL priority subjectId desc organism protein_length source}) {
                    die "No $field for $subjectId" unless $gene->{$field};
                }
                my @pieces = ( a({ -href => $gene->{URL}, -title => $gene->{source} }, $gene->{showName}),
                               $gene->{priority} <= 2 ? b($gene->{desc}) : $gene->{desc},
                               "from",
                               i($gene->{organism}) );
                push @pieces, &simstring(length($seq), $gene->{protein_length},
                                         $queryStart,$queryEnd,$subjectStart,$subjectEnd,
                                         $percIdentity,$eVal,$bitscore,
                                         $def, $gene->{showName}, $seq, $gene->{subjectId})
                    if $gene->{subjectId} eq $genes[0]{subjectId};
                if ($gene->{pmIds}) {
                    my @pmIds = @{ $gene->{pmIds} };
                    my %seen = ();
                    @pmIds = grep { my $keep = !exists $seen{$_}; $seen{$_} = 1; $keep; } @pmIds;
                    my $note = @pmIds > 1 ? scalar(@pmIds) . " papers" : "paper";
                    push @pieces, "(see " .
                        a({-href => "http://www.ncbi.nlm.nih.gov/pubmed/" . join(",",@pmIds) }, $note)
                        . ")";
                }
                push @headers, join(" ", @pieces);
                push @content, $gene->{comment} if $gene->{comment};
                foreach my $paper (@{ $gene->{papers} }) {
                    my @pieces = (); # what to say about this paper
                    my $snippets = [];
                    $snippets = $dbh->selectall_arrayref(
                        "SELECT DISTINCT * from Snippet WHERE geneId = ? AND pmcId = ? AND pmId = ?",
                        { Slice => {} },
                        $gene->{subjectId}, $paper->{pmcId}, $paper->{pmId})
                        if $paper->{pmcId} || $paper->{pmId};
                    my $paperId = join(":::", $paper->{pmId}, $paper->{pmcId}, $paper->{doi});
                    my $nSkip = 0; # number of duplicate snippets
                    foreach my $snippet (@$snippets) {
                        my $text = $snippet->{snippet};
                        my $term = $snippet->{queryTerm};
                        if (exists $paperSeen{$paperId}{$term}) {
                            $nSkip++;
                        } else {
                            $text =~ s!($term)!<B><span style="color: darkred;">$1</span></B>!gi;
                            push @pieces, "...$text...";
                        }
                    }
                    # ignore this paper if all snippets were duplicate terms
                    next if $nSkip == scalar(@$snippets) && $nSkip > 0;
                    foreach my $snippet (@$snippets) {
                        my $term = $snippet->{queryTerm};
                        $paperSeen{$paperId}{$term} = 1;
                    }

                    my $paper_url = undef;
                    my $pubmed_url = "http://www.ncbi.nlm.nih.gov/pubmed/" . $paper->{pmId};
                    if ($paper->{pmcId}) {
                        $paper_url = "http://www.ncbi.nlm.nih.gov/pmc/articles/" . $paper->{pmcId};
                    } elsif ($paper->{pmid}) {
                        $paper_url = $pubmed_url;
                    } elsif ($paper->{doi}) {
                        $paper_url = "http://doi.org/" . $paper->{doi};
                    }
                    my $title = $paper->{title};
                    $title = a({-href => $paper_url}, $title) if defined $paper_url;
                    my $authorShort = $paper->{authors};
                    $authorShort =~ s/ .*//;
                    my $extra = "";
                    $extra = a({-href => $pubmed_url}, "(PubMed)")
                        if !$paper->{pmcId} && $paper->{pmId};
                    my $paper_header = "$title, " .
                        small( a({-title => $paper->{authors}}, "$authorShort,"),
                               $paper->{journal}, $paper->{year}, $extra);
                    
                    if (@$snippets == 0) {
                        # Skip if printed already for this gene (with no snippet)
                        next if exists $paperSeenNoSnippet{$paperId};
                        $paperSeenNoSnippet{$paperId} = 1;

                        # Explain why there is no snippet
                        my $excuse;
                        my $short;
                        if ($paper->{access} eq "full") {
                            $short = "no snippet";
                            $excuse = "This term was not found in the full text, sorry.";
                        } elsif ($paper->{isOpen} == 1) {
                            if ($paper->{access} eq "abstract") {
                                $short = "no snippet";
                                $excuse = "This paper is open access but PaperBLAST only searched the the abstract.";
                            } else {
                                $short = "no snippet";
                                $excuse = "This paper is open access but PaperBLAST did not search either the full text or the abstract.";
                            }
                        } elsif ($paper->{journal} eq "") {
                            $short = "secret";
                            $excuse = "PaperBLAST does not have access to this paper, sorry";
                        } else {
                            $short = "secret";
                            $excuse = "$paper->{journal} is not open access, sorry";
                        }
                        if ($excuse) {
                            my $url = "#secret" if $short eq "secret";
                            my $href = $url ? a({-href => $url, -title => $excuse}, $short)
                                : a({-title => $excuse}, $short);
                            $paper_header .= " " . small("(" . $href . ")"); 
                        }
                    }
                    my $pieces = join(qq{<LI style="list-style-type: none;">}, @pieces);
                    $pieces = qq{<UL style="margin-top: 0em; margin-bottom: 0em;"><LI style="list-style-type: none;">}
                    . $pieces . "</UL>"
                        if $pieces;
                    push @content, $paper_header . $pieces;
                }
            }
            my $content = join("<LI>", @content);
            $content = qq{<UL style="margin-top: 0em; margin-bottom: 0em;"><LI>} . $content . "</UL>"
                if $content;
            print p({-style => "margin-top: 1em; margin-bottom: 0em;"},
                    join("<BR>", @headers) . $content) . "\n";
        }
    }

    print qq{<script src="http://fit.genomics.lbl.gov/d3js/d3.min.js"></script>
             <script src="http://fit.genomics.lbl.gov/images/fitblast.js"></script>
             <H3><A title="Fitness BLAST searches for similarity to bacterial proteins that have mutant phenotypes" HREF="http://fit.genomics.lbl.gov/" NAME="#fitness">Fitness Blast</A></H3>
             <P><DIV ID="fitblast_short"></DIV></P>
             <script>
             var server_root = "http://fit.genomics.lbl.gov/";
             var seq = "$seq";
             fitblast_load_short("fitblast_short", server_root, seq);
             </script>
    };

    my @pieces = $seq =~ /.{1,60}/g;
    print
        h3("Query Sequence"),
        p({-style => "font-family: monospace;"}, small(join(br(), ">$def", @pieces))),
        h3(a({-href => "litSearch.cgi"}, "New Search")),
        $documentation,
        end_html;
}

sub fail($) {
    my ($notice) = @_;
    print
        p($notice),
        p(a({-href => "litSearch.cgi"}, "New search")),
        end_html;
    exit(0);
}

sub simstring($$$$$$$$$) {
    my ($qLen, $sLen, $queryStart,$queryEnd,$subjectStart,$subjectEnd,$percIdentity,$eVal,$bitscore,
        $def1, $def2, $seq1, $acc2) = @_;
    $percIdentity = sprintf("%.0f", $percIdentity);
    # the minimum of coverage either way
    my $cov = ($queryEnd-$queryStart+1) / ($qLen > $sLen ? $qLen : $sLen);
    my $percentCov = sprintf("%.0f", 100 * $cov);
    my $title ="$queryStart:$queryEnd/$qLen of query is similar to $subjectStart:$subjectEnd/$sLen of hit (E = $eVal, $bitscore bits)";
    return "(" .
        a({ -title => $title,
            -href => "showAlign.cgi?def1=$def1&def2=$def2&seq1=$seq1&acc2=$acc2" },
          "$percIdentity% identity, $percentCov% coverage")
        . ")";
}

# The entry will include:
# showName, URL, priority (for choosing what to show first), subjectId, desc, organism, protein_length, source,
# and other entries that depend on the type -- either papers for a list of GenePaper/PaperAccess items,
# or pmIds (a list of pubmed identifiers)
sub SubjectToGene($) {
    my ($subjectId) = @_;
    if ($subjectId =~ m/^gnl\|ECOLI\|(.*)$/) {
        my $protein_id = $1;
        my $gene = $dbh->selectrow_hashref("SELECT * FROM EcoCyc WHERE protein_id = ?", {}, $protein_id);
        die "Unrecognized subject $subjectId" unless defined $gene;
        $gene->{subjectId} = $subjectId;
        $gene->{protein_id} = $protein_id;
        $gene->{source} = "EcoCyc";
        $gene->{URL} = "http://ecocyc.org/gene?orgid=ECOLI&id=$protein_id";
        my @ids = ( $gene->{protein_name}, $gene->{bnumber} );
        @ids = grep { $_ ne "" } @ids;
        $gene->{showName} = join(" / ", @ids) || $protein_id;
        $gene->{priority} = 1;
        $gene->{organism} = "Escherichia coli K-12";
        $gene->{pmIds} = $dbh->selectcol_arrayref(
            "SELECT pmId FROM EcoCycToPubMed WHERE protein_id = ?",
            {}, $protein_id);
        return $gene;
    } else {
        my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE geneId = ?", {}, $subjectId);
        if (defined $gene) {
            $gene->{subjectId} = $subjectId;
            $gene->{priority} = 3;
            my $papers = $dbh->selectall_arrayref(qq{ SELECT DISTINCT * FROM GenePaper
                                                          LEFT JOIN PaperAccess USING (pmcId,pmId)
                                                          WHERE geneId = ? ORDER BY year DESC },
                                                  { Slice => {} }, $subjectId);
            
            $gene->{papers} = $papers;
            my @terms = map { $_->{queryTerm} } @$papers;
            my %terms = map { $_ => 1 } @terms;
            @terms = sort keys %terms;
            $gene->{showName} = join(", ", @terms);
            if ($subjectId =~ m/^VIMSS(\d+)$/) {
                my $locusId = $1;
                $gene->{locusId} = $locusId;
                $gene->{source} = "MicrobesOnline";
                $gene->{URL} = "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$locusId";
                if ($gene->{desc} eq "") {
                    # Fetch the gene description from MicrobesOnline if not in the litsearch db
                    my ($desc) = $mo_dbh->selectrow_array(
                        qq{SELECT description FROM Locus JOIN Description USING (locusId,version)
                           WHERE locusId = ? AND priority=1 },
                        {}, $locusId);
                    $gene->{desc} = $desc || "no description";
                }
            } elsif ($subjectId =~ m/^[A-Z]+_[0-9]+[.]\d+$/) { # refseq
                $gene->{URL} = "http://www.ncbi.nlm.nih.gov/protein/$subjectId";
                $gene->{source} = "RefSeq";
            } elsif ($subjectId =~ m/^[A-Z][A-Z0-9]+$/) { # SwissProt/TREMBL
                $gene->{URL} = "http://www.uniprot.org/uniprot/$subjectId";
                $gene->{source} = "SwissProt/TReMBL";
            } else {
                die "Cannot build a URL for subject $subjectId";
            }
            return $gene;
        } else { # UniProt
            my $gene = $dbh->selectrow_hashref("SELECT * FROM UniProt WHERE acc = ?", {}, $subjectId);
            die "Unrecognized subject $subjectId" unless defined $gene;
            $gene->{subjectId} = $subjectId;
            $gene->{priority} = 2;
            $gene->{source} = "SwissProt";
            $gene->{URL} = "http://www.uniprot.org/uniprot/$subjectId";
            $gene->{showName} = $gene->{acc};
            my @comments = split /_:::_/, $gene->{comment};
            @comments = map { s/[;. ]+$//; $_; } @comments;
            @comments = grep m/^FUNCTION|COFACTOR|CATALYTIC|ENZYME|DISRUPTION/, @comments;
            my $comment = join("<LI>\n", @comments);
            my @pmIds = $comment =~ m!{ECO:0000269[|]PubMed:(\d+)}!;
            $comment =~ s!{ECO:[A-Za-z0-9_:,.| -]+}!!g;
            $gene->{comment} = $comment;
            $gene->{pmIds} = \@pmIds;
            return $gene;
        }
    }
}
