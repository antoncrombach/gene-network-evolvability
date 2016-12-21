//
// Implementation of a facade object, model.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "model.hh"
#include "config.hh"
#include "population.hh"
#include "environment.hh"
#include "observer_manager.hh"
#include "logger.hh"

long blowhole::Model::end_time_ = 0;


blowhole::Model::Model( Fluke *f ) 
    : time_( 0 ), factory_( &( f->configuration() ) ), fluke_( f ),
      poppy_( 0 ), cache_poppy_( 0 ), environ_( 0 ), observers_( 0 ) {
}

blowhole::Model::~Model() {
    if( poppy_ != 0 )
        delete poppy_;
    if( cache_poppy_ != 0 )
        delete cache_poppy_;
    if( observers_ != 0 )
        delete observers_;
    if( environ_ != 0 )
        delete environ_;
}

void
blowhole::Model::build() {
    // build model
#ifdef DEBUG
    cout << "! building environment" << endl;
#endif
    environ_ = factory_.environment();
    environ_->model( this );
#ifdef DEBUG
    cout << "! building population" << endl;
#endif
    cache_poppy_ = factory_.population();
    cache_poppy_->model( this );
    poppy_ = new Population( *cache_poppy_ );
    poppy_->model( this );
    // and how to gather data from it...
#ifdef DEBUG
    cout << "! building observermanager" << endl;
#endif
    observers_ = new ObserverManager();
}

void
blowhole::Model::initialise() {
    // initialise them
    end_time_ = fluke_->configuration().optionAsLong( "end_time" );
    poppy_->initialise();
    environ_->initialise( 
        fluke_->configuration().optionAsInt( "environment_seed" ) );
    observe();
}

void
blowhole::Model::rebuild() {
    // and not to forget
    time_ = 0;
#ifdef DEBUG
    cout << "! building environment" << endl;
#endif 
    if( environ_ != 0 ) {
        delete environ_;
        environ_ = factory_.environment();
        environ_->model( this );
    }
#ifdef DEBUG
    cout << "! reusing population" << endl;
#endif
    if( poppy_ != 0 ) {
        delete poppy_;
        poppy_ = new Population( *cache_poppy_ );
        poppy_->model( this );
    }
#ifdef DEBUG
    cout << "! building observermanager" << endl;
#endif
    if( observers_ != 0 ) {
        delete observers_;
        observers_ = new ObserverManager();
    }
}

void 
blowhole::Model::step() {
    observers_->notifyAll();
    environ_->fluctuate( time_ );
#ifdef DEBUG
std::cout << "!!! time = " << time_ << std::endl;
#endif
    poppy_->step();
    ++time_;
}

bool
blowhole::Model::hasEnded() { 
    return time_ == end_time_ or !poppy_->hasEveryAgentType(); 
}    

void
blowhole::Model::finish() {
    environ_->finish();
    // not needed anymore after change of ancestor logging
    /*poppy_->finish();*/
    unobserve();
}

void
blowhole::Model::observe() {
    // do this in factory?
    Config &aux( fluke_->configuration() );

    if( aux.hasOption( "log_genes_csv" ) or 
        aux.hasOption( "log_active_genes_csv" ) or 
        aux.hasOption( "log_bindingsites_csv" ) or
        aux.hasOption( "log_indegree_csv" ) or 
        aux.hasOption( "log_outdegree_csv" ) or 
        aux.hasOption( "log_ref_genes_csv" ) or 
        aux.hasOption( "log_retroposons_csv" ) ) {
        LogGenomeStats *lgs = new LogGenomeStats( 
            aux.optionAsLong( "log_period_stats" ) );
        observers_->subscribe( poppy_, lgs );
        // and which subobserver do we need?
        if( aux.hasOption( "log_genes_csv" ) ) {
            observers_->subscribe( lgs,
                new LogCsvGenes( aux.optionAsString( "log_genes_csv" ) ) );
        }
        if( aux.hasOption( "log_active_genes_csv" ) ) {
            observers_->subscribe( lgs,
                new LogCsvActiveGenes( aux.optionAsString(
                "log_active_genes_csv" ) ) );
        }
        if( aux.hasOption( "log_bindingsites_csv" ) ) {
            observers_->subscribe( lgs,
                new LogCsvBindingSites( aux.optionAsString(
                "log_bindingsites_csv" ) ) );
        }
        if( aux.hasOption( "log_indegree_csv" ) ) {
            observers_->subscribe( lgs,
                new LogCsvInDegree( aux.optionAsString( "log_indegree_csv" )));
        }
        if( aux.hasOption( "log_outdegree_csv" ) ) {
            observers_->subscribe( lgs,
                new LogCsvOutDegree( aux.optionAsString( 
                    "log_outdegree_csv" ) ) );
        }
        if( aux.hasOption( "log_ref_genes_csv" ) ) {
            observers_->subscribe( lgs,
                new LogCsvRefGenes( aux.optionAsString( 
                    "log_ref_genes_csv" ) ) );
        }
        if( aux.hasOption( "log_retroposons_csv" ) ) {
            observers_->subscribe( lgs,
                new LogCsvRetroposons( aux.optionAsString( 
                    "log_retroposons_csv" ) ) );
        }
    }
    if( aux.hasOption( "log_grid_csv" ) ) {
        observers_->subscribe( poppy_,
            new LogCsvGrid( aux.optionAsString( "log_grid_csv" ),
            aux.optionAsLong( "log_period" ) ) );
    }    
   if( aux.hasOption( "log_genomes_xml" ) ) {
        observers_->subscribe( poppy_, 
            new LogXmlGenomes( aux.optionAsString( "log_genomes_xml" ),
            aux.optionAsLong("log_period_xml" ) ) );
    }
    if( aux.hasOption( "log_distances_csv" ) ) {
        observers_->subscribe( poppy_, 
            new LogCsvDistances( aux.optionAsString( "log_distances_csv" ),
            aux.optionAsLong( "log_period" ) ) );
    }
/*    if( aux.hasOption( "log_perturbations_csv" ) ) {
        observers_->subscribe( poppy_, 
            new LogCsvPerturbations( 
                aux.optionAsString( "log_perturbations_csv" ),
            aux.optionAsLong( "log_period" ) ) );
    }*/
    if( aux.hasOption( "log_scores_csv" ) ) {
        observers_->subscribe( poppy_, 
            new LogCsvScores( aux.optionAsString( "log_scores_csv" ),
            aux.optionAsLong( "log_period" ) ) );
    }
    if( aux.hasOption( "log_population_csv" ) ) {
        observers_->subscribe( poppy_, 
            new LogPopulationSize( aux.optionAsString( "log_population_csv" ),
            aux.optionAsLong( "log_period" ) ) );
    }
    if( aux.hasOption( "log_population_scores_csv" ) ) {
        observers_->subscribe( poppy_, 
            new LogCsvPopulationDistances( 
                aux.optionAsString( "log_population_scores_csv" ),
            aux.optionAsLong( "log_period" ) ) );
    }
    // asynchronous observers
    if( aux.hasOption( "log_environ_csv" ) ) {
        observers_->subscribe( environ_, 
            new LogCsvEnvironment( aux.optionAsString( "log_environ_csv" )));
        environ_->notify();
    }
    // hack!
    if( aux.hasOption( "log_genomes_env_xml" ) ) {
        poppy_->attach2(  
            new LogXmlEnvGenomes( 
                aux.optionAsString( "log_genomes_env_xml" ) ) );
    }
    // more hacking!
    if( aux.hasOption( "log_mutations_csv" ) ) {
        poppy_->attach3(
            new LogCsvMutations( aux.optionAsString( "log_mutations_csv" )));
    }
    // and even more hacks!!`
    if( aux.hasOption( "log_ancestors_csv" ) ) {
        poppy_->attach1(  
            new LogCsvAncestors( aux.optionAsString( "log_ancestors_csv" )));
        // let the logger know its max size
        LogCsvAncestors::maxForestSize( 
            10 * poppy_->grid().shape()[ 0 ] * poppy_->grid().shape()[ 1 ] );
        poppy_->initAncestorTracing();
    } else if( aux.hasOption( "log_agent_trace_xml" ) ) {
        poppy_->attach1(
            new LogXmlAgentTrace( aux.optionAsString( "log_agent_trace_xml" ),
            aux.optionAsString( "agent_trace_source_csv" ) ) );
        poppy_->initAncestorTracing();
    }
}
    
void
blowhole::Model::unobserve() {
    environ_->detachAll();
    observers_->closeAll();
    poppy_->closeAll();
}

