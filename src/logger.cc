//
// Implementation of a few observers that dump their data in files.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "logger.hh"

#include "environment.hh"
#include "population.hh"
#include "genreg_agent.hh"
#include "delta_agent.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "retroposon.hh"
#include "gene.hh"

#include "duo_agent.hh"

//
// Counting double strand breaks and other mutations
//
blowhole::LogCsvMutations::LogCsvMutations( std::string fname ) 
    : AsyncLogObserver( /*s*/ ), time_( -1 ),
     cps_( 4, 0 ), rms_( 4, 0 ), thrs_( 4, 0 ), iacps_( 4, 0 ), 
     iarms_( 4, 0 ), ians_( 4, 0 ), iats_( 4, 0 ), iaws_( 4, 0 ) {
    openLog( fname );
    *log_ << "# WARNING Only handles one module\n";
    *log_ << "# col 1: time\n";
    *log_ << "# col 2-5: pos, neg, neu (cp)\n";
    *log_ << "# col 6-9: pos, neg, neu (rm)\n";
    *log_ << "# col 10-13: pos, neg, neu (thr)\n";
    *log_ << "# col 14-17: pos, neg, neu (ia cp)\n";
    *log_ << "# col 18-21: pos, neg, neu (ia rm)\n";
    *log_ << "# col 22-25: pos, neg, neu (ia nw)\n";
    *log_ << "# col 26-29: pos, neg, neu (ia tag)\n";
    *log_ << "# col 30-34: pos, neg, neu (ia weight)\n";
#ifdef DOUBLESTRANDBREAKS
    *log_ << "# col 35-39: pos, neg, neu (dsb) TODO\n";
#endif
}

void
blowhole::LogCsvMutations::write() {
    if( std::accumulate( cps_.begin(), cps_.end(), 0u ) > 0 or
        std::accumulate( rms_.begin(), rms_.end(), 0u ) > 0 or
        std::accumulate( thrs_.begin(), thrs_.end(), 0u ) > 0 or 
        std::accumulate( iacps_.begin(), iacps_.end(), 0u ) > 0 or
        std::accumulate( iarms_.begin(), iarms_.end(), 0u ) > 0 or
        std::accumulate( ians_.begin(), ians_.end(), 0u ) > 0 or
        std::accumulate( iats_.begin(), iats_.end(), 0u ) > 0 or
        std::accumulate( iaws_.begin(), iaws_.end(), 0u ) > 0 ) {
        // new timestep, write old stuff if nonzero
        *log_ << time_ << "\t\t";
        std::copy( cps_.begin(), cps_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( rms_.begin(), rms_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( thrs_.begin(), thrs_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( iacps_.begin(), iacps_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( iarms_.begin(), iarms_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( ians_.begin(), ians_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( iats_.begin(), iats_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\t";
        std::copy( iaws_.begin(), iaws_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << "\n";
        
        // do not want accumulative: reset to zero
        std::fill( cps_.begin(), cps_.end(), 0 );
        std::fill( rms_.begin(), rms_.end(), 0 );
        std::fill( thrs_.begin(), thrs_.end(), 0 );
        std::fill( iacps_.begin(), iacps_.end(), 0 );
        std::fill( iarms_.begin(), iarms_.end(), 0 );
        std::fill( ians_.begin(), ians_.end(), 0 );
        std::fill( iats_.begin(), iats_.end(), 0 );
        std::fill( iaws_.begin(), iaws_.end(), 0 );
    }
#ifdef DEBUG
    log_->flush();
#endif
}    

void
blowhole::LogCsvMutations::update( Subject *s ) {
    // We are receiving two agents..
    DuoAgent *da = dynamic_cast< DuoAgent* >( s );
    GenRegAgent *m1 = dynamic_cast< GenRegAgent* >( da->first );
    GenRegAgent *m2 = dynamic_cast< GenRegAgent* >( da->second );
    if( m1 and m2 ) {
        // get the time
        if( m1->myTag().time > time_ ) {
            write();
            time_ = m1->myTag().time;
        }
        
        understand( *m1 );
        understand( *m2 );
   }
}

void
blowhole::LogCsvMutations::understand( const GenRegAgent &n ) {
    std::vector< uint > mm( n.nrMutations() );
    // gene dist differences
    int dd = n.deltaDistance();
    int ds = n.deltaSize();

    if( dd < 0 ) {
        // assume something positive happened
        cps_[ POS ] += mm[ Chromosome::CP_G ];
        rms_[ POS ] += mm[ Chromosome::RM_G ];
        thrs_[ POS ] += mm[ Chromosome::THR ];
        iacps_[ POS ] += mm[ Chromosome::CP_IA ];
        iarms_[ POS ] += mm[ Chromosome::RM_IA ];
        ians_[ POS ] += mm[ Chromosome::NW_IA ];
        iats_[ POS ] += mm[ Chromosome::T_IA ];
        iaws_[ POS ] += mm[ Chromosome::W_IA ];
    } else if( dd > 0 ) {
        // negative
        cps_[ NEG ] += mm[ Chromosome::CP_G ];
        rms_[ NEG ] += mm[ Chromosome::RM_G ];
        thrs_[ NEG ] += mm[ Chromosome::THR ];
        iacps_[ NEG ] += mm[ Chromosome::CP_IA ];
        iarms_[ NEG ] += mm[ Chromosome::RM_IA ];
        ians_[ NEG ] += mm[ Chromosome::NW_IA ];
        iats_[ NEG ] += mm[ Chromosome::T_IA ];
        iaws_[ NEG ] += mm[ Chromosome::W_IA ];
    } else {
        if( ds > 0 ) {
            // neutral
            cps_[ NEU ] += mm[ Chromosome::CP_G ];
            iacps_[ NEU ] += mm[ Chromosome::CP_IA ];
            ians_[ NEU ] += mm[ Chromosome::NW_IA ];
        } else if( ds < 0 ) {
            rms_[ NEU ] += mm[ Chromosome::RM_G ];
            iarms_[ NEU ] += mm[ Chromosome::RM_IA ];
        } else {
            thrs_[ NEU ] += mm[ Chromosome::THR ];
            iats_[ NEU ] += mm[ Chromosome::T_IA ];
            iaws_[ NEU ] += mm[ Chromosome::W_IA ];
        }
    }
}

//
// Simple csv observer for the distances
//
blowhole::LogCsvDistances::LogCsvDistances( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
blowhole::LogCsvDistances::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the distances
    std::vector< double > distances;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const GenRegAgent *na = dynamic_cast< GenRegAgent* >( i->first );
        distances.push_back( na->distance() );
    }
    
    // calculate min, median, mean and variance
    if( !distances.empty() ) {
        std::stable_sort( distances.begin(), distances.end() ); 
        double min = distances.front();
        double avg = mean( distances.begin(), distances.end() ); 
        
        // median
        double median = 0.0;
        uint n = distances.size();
        if( n % 2 == 0 ) {
            median = ( distances[ n / 2 ] + distances[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = distances[ ( n / 2 ) ];
        }
        
        // std dev
        double sdev = sqrt( 
            variance( distances.begin(), distances.end(), avg ) );
        
        *log_ << min << "\t" << median << "\t" << avg << "\t" << sdev << "\n";
    } else {
        *log_ << "# Empty grid...\n";
    }
//#ifdef DEBUG
    log_->flush();
//#endif
}

void
blowhole::LogCsvDistances::writeHeader() {
    *log_ << "# min, median, mean, variance of genotypic distances\n";
}

/*//
// Simple csv observer for the perturbations
//
blowhole::LogCsvPerturbations::LogCsvPerturbations( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
blowhole::LogCsvPerturbations::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the Perturbations
    std::vector< double > perturbations;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const GenRegAgent *na = dynamic_cast< GenRegAgent* >( i->first );
        perturbations.push_back( na->nrPerturbations() );
    }
    
    // calculate min, median, mean and variance
    if( !perturbations.empty() ) {
        std::stable_sort( perturbations.begin(), perturbations.end() ); 
        double max = perturbations.back();
        double avg = mean( perturbations.begin(), perturbations.end() ); 
        
        // median
        double median = 0.0;
        uint n = perturbations.size();
        if( n % 2 == 0 ) {
            median = ( perturbations[ n / 2 ] + perturbations[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = perturbations[ ( n / 2 ) ];
        }
        
        // std dev
        double sdev = sqrt( 
            variance( perturbations.begin(), perturbations.end(), avg ) );
        
        *log_ << max << "\t" << median << "\t" << avg << "\t" << sdev << "\n";
    } else {
        *log_ << "# Empty grid...\n";
    }
#ifdef DEBUG
    log_->flush();
#endif
}

void
blowhole::LogCsvPerturbations::writeHeader() {
    *log_ << "# max, median, mean, variance of state perturbations\n";
}*/

//
// Simple csv observer for the configurations
//
blowhole::LogGenomeStats::LogGenomeStats( long i ) 
    : LogObserver( i ), sample_( 0 ), nr_repeats_( 0 ), nr_retroposons_( 0 ), 
        nr_genes_(), nr_genes_on_(), nr_bsites_(), 
        in_degree_(), out_degree_() {}

void
blowhole::LogGenomeStats::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the configurations, bsites and genes
    // not completely sure if this is OOP wise, but it will parallelize better
    nr_genes_.clear();
    nr_genes_on_.clear();
    nr_bsites_.clear();
    in_degree_.clear();
    out_degree_.clear();
    nr_ref_genes_.clear();
    std::fill_n( std::back_inserter( nr_ref_genes_ ), 
        DeltaAgent::nrGenes(), 0 );
    
#ifdef DOUBLESTRANDBREAKS
    nr_repeats_ = 0;
    nr_retroposons_ = 0;
#endif

    sample_ = 0;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const GenRegAgent *na = dynamic_cast< GenRegAgent* >( i->first );
        if( na ) {
            std::list< Chromosome* > chr( na->genome()->chromosomes() );
            // should be only one chromosome here...
            /*for( Genome::chr_iter i = chr.begin(); i != chr.end(); ++i ) {
                doCounting( *i );
            }*/
            doCounting( chr.front() );
            ++sample_;
        }
        // if delta agent, get reference state too
        const DeltaAgent *da = dynamic_cast< DeltaAgent* >( i->first );
        if( da ) {
            doRefState( da );
        }
    }
    // and tell the child observers...
    notify();
}

void
blowhole::LogGenomeStats::doCounting( Chromosome *ch ) {
    // count some things..
    std::set< int > tags;
    std::list< ChromosomeElement* > ces( ch->elements() );
    for( Chromosome::ce_iter i = ces.begin(); i != ces.end(); ++i ) {
        if( Chromosome::IsGene()( *i ) ) {
            Gene *aux = dynamic_cast< Gene* >( *i );
            // gene cp number
            nr_genes_[ aux->tag() ] += 1;
            // active genes
            nr_genes_on_[ aux->tag() ] += aux->state();
            // indegree
            in_degree_[ aux->tag() ] += tags.size();
            // outdegree
            for( std::set< int >::iterator i = tags.begin(); i != tags.end();
                ++i ) {
                out_degree_[ *i ] += 1;
            }
            tags.clear();
        } else if( Chromosome::IsInteraction()( *i ) ) {
            int ia = dynamic_cast< Interaction* >( *i )->tag();
            // bsite cp number
            nr_bsites_[ ia ] += 1;
            // unique interactions per gene
            tags.insert( ia );
        }
#ifdef DOUBLESTRANDBREAKS
          else if( Chromosome::IsRepeat()( *i ) ) {
            ++nr_repeats_;
            /*tags.clear();*/
        } else if( Chromosome::IsRetroposon()( *i ) ) {
            ++nr_retroposons_;
        }
#endif
    }
}

void
blowhole::LogGenomeStats::doRefState( const DeltaAgent *da ) {
    std::vector< int > aux;
    da->referenceState( aux );
    std::transform( nr_ref_genes_.begin(), nr_ref_genes_.end(),
        aux.begin(), nr_ref_genes_.begin(), std::plus< int >() );
}

void
blowhole::LogGenomeStats::nrGenes( std::vector< double > &g ) const {
    int max = 1 + boost::prior( nr_genes_.end() )->first;
    std::fill_n( std::back_inserter( g ), max, 0 );
    for( std::map< int, int >::const_iterator i = nr_genes_.begin(); 
        i != nr_genes_.end(); ++i ) {
        g[ i->first ] = i->second;
    }
}

void
blowhole::LogGenomeStats::nrActiveGenes( std::vector< double > &g ) const {
    int max = 1 + boost::prior( nr_genes_on_.end() )->first;
    std::fill_n( std::back_inserter( g ), max, 0 );
    for( std::map< int, int >::const_iterator i = nr_genes_on_.begin(); 
        i != nr_genes_on_.end(); ++i ) {
        g[ i->first ] = i->second;
    }
}

void
blowhole::LogGenomeStats::nrBindingSites( std::vector< double > &g ) const {
    int max = 1 + boost::prior( nr_bsites_.end() )->first;
    std::fill_n( std::back_inserter( g ), max, 0 );
    for( std::map< int, int >::const_iterator i = nr_bsites_.begin(); 
        i != nr_bsites_.end(); ++i ) {
        g[ i->first ] = i->second;
    }
}

void
blowhole::LogGenomeStats::inDegree( std::vector< double > &g ) const {
    int max = 1 + boost::prior( in_degree_.end() )->first;
    std::fill_n( std::back_inserter( g ), max, 0 );
    for( std::map< int, int >::const_iterator i = in_degree_.begin(); 
        i != in_degree_.end(); ++i ) {
        g[ i->first ] = i->second;
    }
}

void
blowhole::LogGenomeStats::outDegree( std::vector< double > &g ) const {
    int max = 1 + boost::prior( out_degree_.end() )->first;
    std::fill_n( std::back_inserter( g ), max, 0 );
    for( std::map< int, int >::const_iterator i = out_degree_.begin(); 
        i != out_degree_.end(); ++i ) {
        g[ i->first ] = i->second;
    }
}

void
blowhole::LogGenomeStats::nrRefGenes( std::vector< double > &g ) const {
    std::copy( nr_ref_genes_.begin(), nr_ref_genes_.end(), 
        std::back_inserter( g ) );
}

int
blowhole::LogGenomeStats::nrRepeats() const {
    return nr_repeats_;
}

int
blowhole::LogGenomeStats::nrRetroposons() const {
    return nr_retroposons_;
}

int
blowhole::LogGenomeStats::sample() const {
    return sample_;
}

//
// Child observer of GenomeStats
//
blowhole::LogCsvGenes::LogCsvGenes( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# genes\n";
}

void
blowhole::LogCsvGenes::update( Subject *s ) {
    // write averages of each gene cp number
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    double cc = lgs->sample();
    std::vector< double > genes;
    lgs->nrGenes( genes );
    std::transform( genes.begin(), genes.end(), genes.begin(), 
        std::bind2nd( std::divides< double >(), cc ) );
    std::copy( genes.begin(), genes.end(), 
        std::ostream_iterator< double >( *log_, "\t" ) );
    *log_ << "\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 2nd child observer of GenomeStats
//
blowhole::LogCsvBindingSites::LogCsvBindingSites( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# binding sites\n";
}

void
blowhole::LogCsvBindingSites::update( Subject *s ) {
    // write averages of each bsite cp number
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    double cc = lgs->sample();
    std::vector< double > bsites;
    lgs->nrBindingSites( bsites );
    std::transform( bsites.begin(), bsites.end(), bsites.begin(), 
        std::bind2nd( std::divides< double >(), cc ) );
    std::copy( bsites.begin(), bsites.end(), 
        std::ostream_iterator< double >( *log_, "\t" ) );
    *log_ << "\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 3rd child observer of GenomeStats
//
blowhole::LogCsvActiveGenes::LogCsvActiveGenes( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# fraction of genes on\n";
}

void
blowhole::LogCsvActiveGenes::update( Subject *s ) {
    // write averages of each active gene cp number
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    double cc = lgs->sample();
    std::vector< double > genes, genes_on;
    lgs->nrGenes( genes );
    lgs->nrActiveGenes( genes_on );
    std::transform( genes.begin(), genes.end(), genes.begin(), 
        std::bind2nd( std::divides< double >(), cc ) );
    std::transform( genes_on.begin(), genes_on.end(), genes_on.begin(), 
        std::bind2nd( std::divides< double >(), cc ) );
    std::transform( genes_on.begin(), genes_on.end(), genes.begin(),
        genes.begin(), std::divides< double >() );
    std::copy( genes.begin(), genes.end(), 
        std::ostream_iterator< double >( *log_, "\t" ) );
    *log_ << "\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 4th child observer of GenomeStats
//
blowhole::LogCsvInDegree::LogCsvInDegree( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# nr of unique bsites per gene\n";
}

void
blowhole::LogCsvInDegree::update( Subject *s ) {
    // write averages of each genes in degree (lumped by tag)
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    double cc = lgs->sample();
    std::vector< double > targets;
    lgs->inDegree( targets );
    std::transform( targets.begin(), targets.end(), targets.begin(), 
        std::bind2nd( std::divides< double >(), cc ) );
    std::copy( targets.begin(), targets.end(), 
        std::ostream_iterator< double >( *log_, "\t" ) );
    *log_ << "\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 5th child observer of GenomeStats
//
blowhole::LogCsvOutDegree::LogCsvOutDegree( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# nr of unique targets per gene\n";
}

void
blowhole::LogCsvOutDegree::update( Subject *s ) {
    // write averages of each genes out degree (lumped by tag)
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    double cc = lgs->sample();
    std::vector< double > targets;
    lgs->outDegree( targets );
    std::transform( targets.begin(), targets.end(), targets.begin(), 
        std::bind2nd( std::divides< double >(), cc ) );
    std::copy( targets.begin(), targets.end(), 
        std::ostream_iterator< double >( *log_, "\t" ) );
    *log_ << "\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 6th child observer of GenomeStats
//
blowhole::LogCsvRefGenes::LogCsvRefGenes( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# nr of unique targets per gene\n";
}

void
blowhole::LogCsvRefGenes::update( Subject *s ) {
    // write averages of each gene's reference state 
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    double cc = lgs->sample();
    std::vector< double > targets;
    lgs->nrRefGenes( targets );
    std::transform( targets.begin(), targets.end(), targets.begin(), 
        std::bind2nd( std::divides< double >(), cc ) );
    std::copy( targets.begin(), targets.end(), 
        std::ostream_iterator< double >( *log_, "\t" ) );
    *log_ << "\n";
#ifdef DEBUG
    log_->flush();
#endif
}

//
// 7th child observer of GenomeStats
//
blowhole::LogCsvRetroposons::LogCsvRetroposons( 
        std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    *log_ << "# retroposons repeats\n";
}

void
blowhole::LogCsvRetroposons::update( Subject *s ) {
    // write averages of each gene's reference state 
    LogGenomeStats *lgs = dynamic_cast< LogGenomeStats* >( s );
    double cc = lgs->sample();
    *log_ << lgs->nrRetroposons() / cc << "\t" << lgs->nrRepeats() / cc << "\n";
#ifdef DEBUG
    log_->flush();
#endif
}


//
// Simple csv observer for the scores
//
blowhole::LogCsvScores::LogCsvScores( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
blowhole::LogCsvScores::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the scores
    std::vector< double > scores;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        scores.push_back( i->first->score() );
    }
    
    // calculate max, median, mean and variance
    if( !scores.empty() ) {
        std::stable_sort( scores.begin(), scores.end() ); 
        double max = scores.back();
        double avg = mean( scores.begin(), scores.end() ); 
        
        // median
        double median = 0.0;
        uint n = scores.size();
        if( n % 2 == 0 ) {
            median = ( scores[ n / 2 ] + scores[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = scores[ ( n / 2 ) ];
        }
        
        // std dev
        double sdev = sqrt( 
            variance( scores.begin(), scores.end(), avg ) );
       
        *log_ << max << "\t" << median << "\t" << avg << "\t" << sdev << "\n";
    } else {
        *log_ << "# Empty grid...\n";
    }
    log_->flush();
}

void
blowhole::LogCsvScores::writeHeader() {
    *log_ << "# max, median, mean, variance of raw fitness scores\n";
}


//
// Simple csv observer for the environment
//
blowhole::LogCsvEnvironment::LogCsvEnvironment( 
        std::string fname ) : AsyncLogObserver( /*s*/ ) {
    openLog( fname );
    *log_ << "# time, new environment\n";
}

void
blowhole::LogCsvEnvironment::update( Subject *s ) {
    Environment *env = dynamic_cast< ConstantEnvironment * >( s );
    if( !env ) env = dynamic_cast< PoissonEnvironment * >( s );
    if( !env ) env = dynamic_cast< PeriodicEnvironment * >( s );
    // log change from which to which state
    *log_ << env->model()->now() << "\t" << env->attractor() << "\n";
    log_->flush();
}


//
// Simple csv observer for ancestor tracing
//
/*
blowhole::LogCsvAncestors::LogCsvAncestors( 
    std::string fname ) : AsyncLogObserver() {
    openLog( fname );
    writeHeader();
}

void 
blowhole::LogCsvAncestors::update( Subject *s ) {
    Agent *ag = dynamic_cast< Agent * >( s );
    AgentTag child = ag->myTag();
    AgentTag mother = ag->parentTag();
    // log birth of new agent
    *log_ << child.time << "\t" << child.loc.first << "\t" 
        << child.loc.second << "\t" << child.i << "\t-> ";
    *log_ << mother.time  << "\t" << mother.loc.first << "\t"
        << mother.loc.second << "\t" << mother.i << "\n";
    log_->flush();
}

void 
blowhole::LogCsvAncestors::writeHeader() {
    *log_ << "# child ( birth, x, y ) -> parent ( birth, x, y )\n";
}

void 
blowhole::LogCsvAncestors::initialize( Subject *s ) {}

void
blowhole::LogCsvAncestors::finalize() {}
*/

//
// More complicated csv observer for ancestor tracing
//
// Note: have to make an educated guess for the max allowed size,
// should depend on grid size...
uint blowhole::LogCsvAncestors::max_forest_size_ = 100000;

blowhole::LogCsvAncestors::LogCsvAncestors( 
    std::string fname ) : AsyncLogObserver(), nr_roots_( 0 ), forest_(),
    leafs_(), mothers_() {
    openLog( fname );
    writeHeader();
    initForest();
}

void
blowhole::LogCsvAncestors::initialize( Subject *s ) {
    // init agent by putting its id in the leafs and forest
    Agent *ag = dynamic_cast< Agent * >( s );
    AgentTag ch = ag->myTag();
    AgentTag mo;
    
    leafs_.insert( std::make_pair( ch, boost::add_vertex( forest_ ) ) );
    boost::put( boost::vertex_bundle, forest_, leafs_[ ch ], ch );
    // \c mo is null object and present in forest after initForest
    boost::add_edge( leafs_[ ch ], mothers_[ mo ], forest_ );
    ++nr_roots_;
}

void 
blowhole::LogCsvAncestors::update( Subject *s ) {
    Agent *ag = dynamic_cast< Agent * >( s );
    if( nr_roots_ == 1 ) {
        // renew buffer, write old contents
        writeForest();
        renewForest();
    } 
    if( boost::num_vertices( forest_ ) > max_forest_size_ ) {
        // if forest becomes really large...
        writeForest();
        renewForest();
    }
    // what to do?
    if( leafs_.find( ag->myTag() ) != leafs_.end() ) {
        death( ag );
    } else {
        birth( ag );
    }
}

void
blowhole::LogCsvAncestors::birth( Agent *ag ) {
    AgentTag ch = ag->myTag();
    AgentTag mo = ag->parentTag();

    std::map< AgentTag, anc_vertex >::iterator aux = leafs_.find( mo );
    if( aux != leafs_.end() ) {
        // not a mother yet, check grandma too (if there is any...)
        AgentTag g;
        adj_iter vg, end;
        boost::tie( vg, end ) = boost::adjacent_vertices( aux->second,forest_ );
        if( vg != end ) { 
            g = boost::get( boost::vertex_bundle, forest_, *vg );
        }
        if( !g.isNull() ) {
            invadj_iter i, j, k;
            boost::tie( i, j ) = boost::inv_adjacent_vertices( *vg, forest_ );
            k = j;
            while( i != j ) {
                std::map< AgentTag, anc_vertex >::iterator bux =
                    leafs_.find( 
                        boost::get( boost::vertex_bundle, forest_, *i ) );
                if( bux == aux or bux == leafs_.end() ) { 
                    ++i;
                } else {
                    j = i;
                }
            }
            if( j == k ) {
                // no more leafs to be found, grandma shouldn't be a mom
                mothers_.erase( g );
            }
        }
        // add mo to mother and remove it from leafs
        mothers_.insert( *aux );
        leafs_.erase( aux );
    } else if( mothers_.find( mo ) == mothers_.end() ) {
        // only happens if inbetween two births the forest has been written
        AgentTag rt;
        mothers_.insert( std::make_pair( mo, boost::add_vertex( forest_ ) ) );
        boost::put( boost::vertex_bundle, forest_, mothers_[ mo ], mo );
        boost::add_edge( mothers_[ mo ], mothers_[ rt ], forest_ );
    }

    // and add the child
    leafs_.insert( std::make_pair( ch, boost::add_vertex( forest_ ) ) );
    boost::put( boost::vertex_bundle, forest_, leafs_[ ch ], ch );
    boost::add_edge( leafs_[ ch ], mothers_[ mo ], forest_ );
}

void
blowhole::LogCsvAncestors::death( Agent *ag ) {
    AgentTag ch = ag->myTag();
    AgentTag mo = ag->parentTag();

    boost::clear_vertex( leafs_[ ch ], forest_ );
    boost::remove_vertex( leafs_[ ch ], forest_ );
    leafs_.erase( ch );
    // mother is special case ( if not isNull )
    anc_vertex j, i;
    std::map< AgentTag, anc_vertex >::iterator mother( mothers_.find( mo ) );
    if( mother != mothers_.end() ) {
        i = mother->second;
        if( i != root_ ) {
            if( boost::in_degree( i, forest_ ) == 0 ) {
                j = *( boost::adjacent_vertices( i, forest_ ).first );
                boost::clear_vertex( i, forest_ );
                boost::remove_vertex( i, forest_ );
                mothers_.erase( mother );
                i = j;
            } else {
                // other possibility is in_degree == 1
                // if other branch is not a leaf, mother is no more a mother 
                j = *( boost::inv_adjacent_vertices( i, forest_ ).first );
                if( boost::in_degree( j, forest_ ) != 0 ) {
                    mothers_.erase( mother );
                }
            }
        }  
    } else {
        i = root_;
    }

    // further down the lineage
    while( i != root_ and boost::in_degree( i, forest_ ) == 0 ) {
        j = *( boost::adjacent_vertices( i, forest_ ).first );
        boost::clear_vertex( i, forest_ );
        boost::remove_vertex( i, forest_ );
        i = j;
    }
    // did we hit a root?
    if( i == root_ ) {
        --nr_roots_;
    }
}

void
blowhole::LogCsvAncestors::finalize() {
    writeForest();
#ifdef DEBUG
    *log_ << "}";
#endif
}

void 
blowhole::LogCsvAncestors::writeHeader() {
    *log_ << "# child ( birth, x, y ) -> parent ( birth, x, y )\n";
#ifdef DEBUG
    *log_ << "digraph G {\n";
#endif
}

void 
blowhole::LogCsvAncestors::writeForest() {
    // log part of the total ancestor tree, by time;
    std::map< AgentTag, anc_vertex > time_map;
    mapToTime( time_map );
    // according to map, write to log
    for( std::map< AgentTag, anc_vertex >::iterator i = time_map.begin();
        i != time_map.end(); ++i ) {
        AgentTag mother =  boost::get( boost::vertex_bundle, forest_, 
            *( boost::adjacent_vertices( i->second, forest_ ).first ) );
        if( !i->first.isNull() and !mother.isNull() ) {
#ifdef DEBUG
            *log_ << "\"" << i->first << "\"\t-> \"" << mother << "\"\n";
#else
            *log_ << i->first << "\t-> " << mother << "\n";
#endif
        }
    }
#ifdef DEBUG
    for( std::map< AgentTag, anc_vertex >::iterator i = mothers_.begin();
        i != mothers_.end(); ++i ) {
        *log_ << "\"" << i->first << "\"[color=\"blue\"];\n";
    }
    for( std::map< AgentTag, anc_vertex >::iterator i = leafs_.begin();
        i != leafs_.end(); ++i ) {
        *log_ << "\"" << i->first << "\"[shape=\"diamond\"];\n";
    }
#endif

#ifdef DEBUG
    log_->flush();
#endif
}

void
blowhole::LogCsvAncestors::renewForest() {
    // pre: one root?
    // empty graph
    forest_.clear();
    mothers_.clear();
    // add null object
    AgentTag aux;
    root_ = boost::add_vertex( forest_ );
    boost::put( boost::vertex_bundle, forest_, root_, aux );
    mothers_.insert( std::make_pair( aux, root_ ) );
    // and leafs at the start..
    for( std::map< AgentTag, anc_vertex >::iterator i = leafs_.begin();
        i != leafs_.end(); ++i ) {
        i->second = boost::add_vertex( forest_ );
        boost::put( boost::vertex_bundle, forest_, i->second, i->first );
        boost::add_edge( i->second, root_, forest_ );
    }
    nr_roots_ = leafs_.size();
}

void
blowhole::LogCsvAncestors::initForest() {
    // init forest by adding the ultimate root to the mothers section
    AgentTag aux;
    root_ = boost::add_vertex( forest_ );
    boost::put( boost::vertex_bundle, forest_, root_, aux );
    mothers_.insert( std::make_pair( aux, root_ ) );    
}

void
blowhole::LogCsvAncestors::mapToTime( std::map< AgentTag, anc_vertex > &tm ) {
    anc_birth_map tag_map = boost::get( boost::vertex_bundle, forest_ );
    ver_iter n, m;
    for( boost::tie( n, m ) = boost::vertices( forest_ ); n != m; ++n ) {
        tm.insert( std::make_pair( tag_map[ *n ], *n ) );
    }
}

void
blowhole::LogCsvAncestors::maxForestSize( uint mfs ) 
{ max_forest_size_ = mfs; }

//
// After one run, we can trace back agents to the beginning and then follow
// their development
//
blowhole::LogXmlAgentTrace::LogXmlAgentTrace( 
    std::string dname, std::string tracefile ) : AsyncLogObserver() {
    first_ = true;
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
    trace_ = StreamManager::instance()->openInFileStream( 
        tracefile, std::fstream::in );
    nextTarget();
}

void 
blowhole::LogXmlAgentTrace::initialize( Subject *s ) { 
    update( s );
}

void
blowhole::LogXmlAgentTrace::update( Subject *s ) {
    GenRegAgent *ag = dynamic_cast< GenRegAgent * >( s );
    AgentTag tag = ag->myTag();
    if( tag == target_ ) {
        // devise name
        std::string fname = tag.str();
        // write genome
        openLog( dname_ + "/" + fname + ".xml" );
        writeHeader();
#ifdef DOUBLESTRANDBREAKS
        std::vector< uint > aux( ag->nrMutations() );
        *log_ << "<mutations dsb=\"" << aux[ Chromosome::DSB ]
            << "\" cp_g=\"" << aux[ Chromosome::CP_G ] 
            << "\" rm_g=\"" << aux[ Chromosome::RM_G ] << "\"/>\n";
#endif
        *log_ << *ag;
        writeFooter();
        closeLog();
        // and move on
        nextTarget();
    }
}

void
blowhole::LogXmlAgentTrace::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation fluke_version=\"" << VERSION << "\">\n";
}

void
blowhole::LogXmlAgentTrace::writeFooter() {
    *log_ << "</simulation>\n";
}

void
blowhole::LogXmlAgentTrace::nextTarget() {
    // checks for input / output necessary?
    *trace_ >> target_.time >> target_.loc.first >> target_.loc.second 
            >> target_.i;
}


//
// Another class
// 
blowhole::LogXmlGenomes::LogXmlGenomes( std::string dname, long i ) 
    : LogObserver( i ) {
    // create dir with given name
    // within dir, use some logical name
    // f.i. timestep.dot
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
}

void
blowhole::LogXmlGenomes::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    
    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    *log_ << "<generation time=\"" << pop->generation() << "\">\n"
          << *pop << "</generation>\n";
    closeLog();
}

void
blowhole::LogXmlGenomes::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation fluke_version=\"" << VERSION << "\">\n";
}

void
blowhole::LogXmlGenomes::writeFooter() {
    *log_ << "</simulation>\n";
}

std::string
blowhole::LogXmlGenomes::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".xml";
    return result.str();
}

//
// Another class, yet asynchronous
// 
blowhole::LogXmlEnvGenomes::LogXmlEnvGenomes( std::string dname ) 
    : AsyncLogObserver() {
    // create dir with given name
    // within dir, use some logical name
    // f.i. timestep_xcoordinate_ycoordinate.dot
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
}

void
blowhole::LogXmlEnvGenomes::update( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
   
    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    *log_ << "<generation time=\"" << pop->generation() << "\">\n"
          << *pop << "</generation>\n";
    //writeFooter();
    closeLog();
}

void
blowhole::LogXmlEnvGenomes::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation fluke_version=\"" << VERSION << "\">\n";
}

void
blowhole::LogXmlEnvGenomes::writeFooter() {
    *log_ << "</simulation>\n";
}

std::string
blowhole::LogXmlEnvGenomes::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".xml";
    return result.str();
}


//
// Yet another class
//
blowhole::LogPopulationSize::LogPopulationSize(
    std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
blowhole::LogPopulationSize::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    *log_ << pop->nrAgents() << std::endl;
}

void
blowhole::LogPopulationSize::writeHeader() {
    *log_ << "# population size\n";
}

//
// Simple csv observer for the population scores
//
blowhole::LogCsvPopulationDistances::LogCsvPopulationDistances( 
        std::string fname, long i ) : LogObserver( i ) {
    openLog( fname );
    writeHeader();
}

void
blowhole::LogCsvPopulationDistances::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the scores
    std::vector< double > scores_one;
    std::vector< double > scores_two;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        const GenRegAgent *na = dynamic_cast< GenRegAgent* >( i->first );
        if( na->type() == 1 ) {
            scores_one.push_back( na->distance() );
        } else if( na->type() == 2 ) {
            scores_two.push_back( na->distance() );
        }
    }
    
    // calculate min, median, mean and variance for score_one
    if( !scores_one.empty() ) {
        std::stable_sort( scores_one.begin(), scores_one.end() ); 
        double min = scores_one.front();
        double avg = mean( scores_one.begin(), scores_one.end() ); 
        
        // median
        double median = 0.0;
        uint n = scores_one.size();
        if( n % 2 == 0 ) {
            median = ( scores_one[ n / 2 ] + scores_one[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = scores_one[ ( n / 2 ) ];
        }
        
        // std dev
        double sdev = sqrt( 
            variance( scores_one.begin(), scores_one.end(), avg ) );
        
        *log_ << "1\t" << min << "\t" << median 
              << "\t" << avg << "\t" << sdev << "\t";
    } else {
        *log_ << "1\t0\t0\t0\t0\t";
    }
    // calculate min, median, mean and variance
    if( !scores_two.empty() ) {
        std::stable_sort( scores_two.begin(), scores_two.end() ); 
        double min = scores_two.front();
        double avg = mean( scores_two.begin(), scores_two.end() ); 
        
        // median
        double median = 0.0;
        uint n = scores_two.size();
        if( n % 2 == 0 ) {
            median = ( scores_two[ n / 2 ] + scores_two[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = scores_two[ ( n / 2 ) ];
        }
        
        // std dev
        double sdev = sqrt( 
            variance( scores_two.begin(), scores_two.end(), avg ) );
        
        *log_ << "2\t" << min << "\t" << median 
              << "\t" << avg << "\t" << sdev << "\n";
    } else {
        *log_ << "\n";
    }
    log_->flush();
}

void
blowhole::LogCsvPopulationDistances::writeHeader() {
    *log_ << "# min, median, mean, variance of raw scores per agent type\n";
}

//
// And now a population dumper, at least it dumps a few specific feats
//
blowhole::LogCsvGrid::LogCsvGrid( std::string dname, long i ) 
    : LogObserver( i ) {
    dname_ = dname;
    StreamManager::instance()->openPath( dname_ );
}

void
blowhole::LogCsvGrid::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_grid grid = pop->grid();
    
    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    // and a label of the feat
    *log_ << "gene distance" << "\t";
    // with the dimensions of the matrix to follow
    int n = grid.shape()[ 0 ];
    int m = grid.shape()[ 1 ];
    *log_ << n << "\t" << m << "\n";
    for( int i = 0; i != n; ++i ) {
        for( int j = 0; j != m; ++j ) {
            if( grid[ i ][ j ] != 0 ) {
                *log_ << grid[ i ][ j ]->score() << "\t";
            } else {
                *log_ << -1 << "\t";
            }
        }
        *log_ << "\n";
    }
    closeLog();
}

void
blowhole::LogCsvGrid::writeHeader() {
    *log_ << "# grid size and features in matrix format,";
    *log_ << " separate by newlines\n";
}

std::string
blowhole::LogCsvGrid::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".csv";
    return result.str();
}

