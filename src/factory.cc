//
// Implementation of the model factory
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "factory.hh"
#include "stream_manager.hh"

#include "environment.hh"
#include "population.hh"
#include "simple_agent.hh"
#include "net_agent.hh"
#include "delta_agent.hh"

#include "genome.hh"
#include "chromosome.hh"
#include "gene.hh"
#include "interaction.hh"
#include "retroposon.hh"
#include "repeat.hh"

#include "hop_graph.hh"

#include "scaling.hh"
#include "selection.hh"

#include "net_agent_reader.hh"
#include "net_pop_reader.hh"

XERCES_CPP_NAMESPACE_USE

blowhole::Factory::Factory() {
    gene_tag_ = 0;
}

blowhole::Factory::Factory( Config *gc ) : conf_( gc ) {
    gene_tag_ = 0;
}

blowhole::Factory::~Factory() {}

blowhole::Gene*
blowhole::Factory::gene() {
    return new Gene( tag() );
}

blowhole::Retroposon* 
blowhole::Factory::retroposon() {
    return new Retroposon( tag() );
}

blowhole::Repeat*
blowhole::Factory::repeat() {
    return new Repeat();
}

blowhole::Chromosome* 
blowhole::Factory::chromosome( int type ) {
    // auxilary typedef
    std::list< ChromosomeElement* > *ll = 
        organiseChromosome( chromosomeParts( type ), 
            conf_->optionAsDouble( "organised", type ) );
    
    // scan list and flank all retroposons with repeats
    Chromosome::ce_iter i = ll->begin();
    while( i != ll->end() ) {
        if( Chromosome::IsRetroposon()( *i ) ) {
            ll->insert( i, repeat() );
            ll->insert( boost::next( i ), repeat() );
        }
        ++i;
    }
    
    // init with right values for this type of agent
    Chromosome *result = new Chromosome( 0, ll );
    return result;
}

blowhole::Genome* 
blowhole::Factory::genome( int type ) {
    resetTag();
    // #chromosomes
    int n_ch = conf_->optionAsInt( "chromos", type );

    std::list< Chromosome* > *ll = new std::list< Chromosome* >();
    for( int i = 0; i < n_ch; ++i ) {
        ll->push_back( this->chromosome( type ) );
    }
    return new Genome( ll );   
}

void
blowhole::Factory::network( int type ) {
#ifdef DEBUG
    cout << "! init reference network" << endl;
#endif
    // read reference network
    boost::filesystem::ifstream *file = 
        StreamManager::instance()->openInFileStream( conf_->optionAsString( "network", type ), boost::filesystem::fstream::in );
    HopfieldGraph::readReferenceGraph( *file );
    StreamManager::instance()->closeInFileStream( file );
    HopfieldGraph::maxPropagate( conf_->optionAsInt( "max_propagate", type ) );
    HopfieldGraph::seqPropagate( boost::to_lower_copy( 
        conf_->optionAsString( "seq_propagate", type ) ) == "true" );
    HopfieldGraph::perturbRate( conf_->optionAsDouble( "perturb_rate", type ) );

    // write reference network
    boost::filesystem::ofstream *outfile = StreamManager::instance()->
        openOutFileStream( std::string( "refnet.dot" ),
        boost::filesystem::fstream::out );
    HopfieldGraph::writeReferenceGraph( *outfile );
    StreamManager::instance()->closeOutFileStream( outfile );
}

blowhole::Agent*
blowhole::Factory::agent( int type ) {
    // build agent according to type    
    Agent *result;
    if( conf_->optionAsString( "agent", type ) == "simple" ) {
        result = new SimpleAgent();
    } else if( conf_->optionAsString( "agent", type ) == "net" ) {
        result = new NetAgent( type, genome( type ) );
        NetAgent::sensor( conf_->optionAsInt( "env_bit", type ) );
        network( type );
    } else if( conf_->optionAsString( "agent", type ) == "delta" ) {
        boost::dynamic_bitset<> aux( 
            conf_->optionAsString( "init_ref_state", type ) );
        result = new DeltaAgent( type, genome( type ), aux );
        DeltaAgent::sensor( conf_->optionAsInt( "env_bit", type ) );
        network( type );
    } else {
        // try to open it as a file
        result = readNetAgent( conf_->optionAsString( "agent", type ), type );
        network( type );
    }
    result->initialise();
    return result;
}

blowhole::Population*
blowhole::Factory::population() {
    Population *result = 0;
    // read some population level parameters
    Population::nrAgentTypes( conf_->optionAsInt( "nr_agent_type" ) );
    Population::placement( conf_->optionAsString( "agent_placement" ) );
    Population::shuffling( conf_->optionAsString( "shuffle" ) == "true" );
    /*
    Population::threshold( conf_->optionAsDouble( "sum_fitness_threshold" ) );
    */
    // and per agent type stuff
    readAgentConfigurations();
    
    uint xx = conf_->optionAsInt( "grid_x" );
    uint yy = conf_->optionAsInt( "grid_y" );
    bool first_pop = conf_->hasOption( "population_one" );
    bool second_pop = conf_->hasOption( "population_two" );
    bool rand_pop = conf_->optionAsString( "population_start" ) == "random"; 
    if( !first_pop and !second_pop ) {
        std::cout << "Creating a population" << std::endl;
        std::vector< Agent* > aux;
        int cux = conf_->optionAsInt( "nr_agent_type" );
        // random start or homogeneous start
        if( rand_pop ) {
            for( int i = 0; i != cux; ++i ) {
                int bux = conf_->optionAsInt( "init_nr_agents", i );
                for( int j = 0; j != bux; ++j ) {
                    aux.push_back( agent( i ) );
                }
            }
        } else {
            for( int i = 0; i != cux; ++i ) {
                int bux = conf_->optionAsInt( "init_nr_agents", i );
                if( bux > 0 ) {
                    aux.push_back( agent( i ) );
                    for( int j = 1; j != bux; ++j ) {
                        aux.push_back( aux.front()->clone() );
                    }
                }
            }
        }
        result = new Population( xx, yy, aux,
            scalingScheme(), selectionScheme() );
    } else if( first_pop and !second_pop ) {
        std::cout << "Reading in one population (1)" << std::endl;
        std::vector< Agent* > aux;
        std::vector< location > loc;
        boost::tie( aux, loc ) = 
            readNetPop( conf_->optionAsString( "population_one" ), 0 );
        // set agent type to 0 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            //( **ii ).type( 0 );
            ( **ii ).initialise();
        } 
        result = new Population( xx, yy, aux, loc,
            scalingScheme(), selectionScheme() );
    } else if( !first_pop and second_pop ) {
        std::cout << "Reading in one population (2)" << std::endl;
        std::vector< Agent* > aux;
        std::vector< location > loc;
        boost::tie( aux, loc ) = 
            readNetPop( conf_->optionAsString( "population_two" ), 0 );
        // set agent type to 0 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            //( **ii ).type( 0 );
            ( **ii ).initialise();
        } 
        result = new Population( xx, yy, aux, loc,
            scalingScheme(), selectionScheme() );
    } else { // both populations
        std::cout << "Reading in two populations" << std::endl;
        std::vector< Agent* > aux;
        std::vector< location > loc;
        boost::tie( aux, loc ) = 
            readNetPop( conf_->optionAsString( "population_one" ), 0 );
        // set agent type to 0 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            ( **ii ).type( 0 );
            ( **ii ).initialise();
        } 
        // 2nd population
        std::vector< Agent* > bux;
        std::vector< location > bloc;
        boost::tie( bux, bloc ) = 
            readNetPop( conf_->optionAsString( "population_two" ), 1 );
        // translation of locations
        uint dx = conf_->optionAsInt( "grid_x" ) / 2;
        for( Population::loc_iter ii = bloc.begin();
            ii != bloc.end(); ++ii ) {
            ii->first += dx;
        }
        // set agent type to 1 and init
        for( std::vector< Agent* >::iterator ii = bux.begin(); 
            ii != bux.end(); ++ii ) {
            ( **ii ).type( 1 );
            ( **ii ).initialise();
        }
        result = new Population( xx, yy, bux, bloc,
            scalingScheme(), selectionScheme() );
    }
    return result;
}

blowhole::Environment*
blowhole::Factory::environment() {
    // FIX using magic numbers for the number of states
    Environment *result;
    if( conf_->optionAsString( "environment" ) == "constant" ) {
        result = new ConstantEnvironment( 1 );
        ConstantEnvironment *aux =
            dynamic_cast< ConstantEnvironment * >( result );
        // reading in string of bits
        std::vector< boost::dynamic_bitset<> > bux;
        readEnvironment( bux );
        aux->attractor( bux[ 0 ] );
    } else if( conf_->optionAsString( "environment" ) == "poisson" ) {
        // get attractors from configuration
        std::vector< boost::dynamic_bitset<> > bux;
        readEnvironment( bux );
        result = new PoissonEnvironment( bux.size() );
        PoissonEnvironment *aux =
            dynamic_cast< PoissonEnvironment * >( result );
        for( uint i = 0; i < bux.size(); ++i ) {
            aux->attractor( i, bux[ i ] );
        }        
        aux->lambda( conf_->optionAsDouble( "lambda" ) );
    } else if( conf_->optionAsString( "environment" ) == "periodic" ) {
        result = new PeriodicEnvironment( 2 );
        PeriodicEnvironment *aux =
            dynamic_cast< PeriodicEnvironment * >( result );
        // reading in string of bits
        std::vector< boost::dynamic_bitset<> > bux;
        readEnvironment( bux );
        // and assign them
        aux->attractor( 0, bux[ 0 ] );
        aux->attractor( 1, bux[ 1 ] );
        aux->period( static_cast< long >( conf_->optionAsDouble( "lambda" ) ) );
        aux->offset( conf_->optionAsLong( "offset" ) );
    } else {
        result = new ConstantEnvironment( 1 );
        ConstantEnvironment *aux =
            dynamic_cast< ConstantEnvironment * >( result );
        // reading in string of bits
        std::vector< boost::dynamic_bitset<> > bux;
        readEnvironment( bux );
        aux->attractor( bux[ 0 ] );
    }
    return result;
}

void
blowhole::Factory::readEnvironment( 
    std::vector< boost::dynamic_bitset<> > &attr ) {
    std::vector< std::string > bux( conf_->optionAsVector( "attractors" ) );
    // reading in string of bits
    for( std::vector< std::string >::iterator i = bux.begin(); i != bux.end();
        ++i ) {
        try {
            // decimal
            boost::dynamic_bitset<> cux( 8 * sizeof( ulong ), 
                boost::lexical_cast< ulong >( *i ) );
            attr.push_back( cux );
        } catch( boost::bad_lexical_cast &bc ) {
            // binary (and remove the '0b' header)
            boost::dynamic_bitset<> cux( i->erase( 0, 2 ) );
            attr.push_back( cux );
        }
    }
}

blowhole::ScalingScheme*
blowhole::Factory::scalingScheme() {
    ScalingScheme *result;
    if( conf_->optionAsString( "scaling_scheme" ) == "none" ) {
        result = new NoScaling();
        Population::threshold( THRESHOLD );
    } else if( conf_->optionAsString( "scaling_scheme" ) == "linear" ) {
        result = new LinearScaling( conf_->optionAsDouble( "base_score" ) );
        Population::threshold( THRESHOLD );
    } else if( conf_->optionAsString( "scaling_scheme" ) == "power" ) {
        result = new PowerScaling( conf_->optionAsDouble( "base_score" ), 
            SELECTION );
        Population::threshold( std::pow( THRESHOLD, SELECTION ) );
    } else if( conf_->optionAsString( "scaling_scheme" ) == "exponential" ) {
        result = new ExponentialScaling( conf_->optionAsDouble( "base_score" ),
            SELECTION );
        Population::threshold( std::exp( THRESHOLD * SELECTION ) );
    } else {
        // default
        result = new NoScaling();
        Population::threshold( THRESHOLD );
    }
    return result;
}

blowhole::SelectionScheme*
blowhole::Factory::selectionScheme() {
    SelectionScheme *result;
    if( conf_->optionAsString( "selection_scheme" ) == "random" ) {
        result = new RandSelection();
    } else if( conf_->optionAsString( "selection_scheme" ) == "probalistic" ) {
        result = new ProbalisticSelection();
    } else {
        // default
        result = new RandSelection();
    }
    return result;
}


blowhole::GenRegAgent*
blowhole::Factory::readNetAgent( std::string fname, int type ) {
    // Part for reading the genome
    try {
        XMLPlatformUtils::Initialize();
    } catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Error during initialization! :\n" << message << "\n";
        XMLString::release( &message );
    }

    SAX2XMLReader *parser = XMLReaderFactory::createXMLReader();
    parser->setFeature( XMLUni::fgSAX2CoreValidation, false );
    parser->setFeature( XMLUni::fgSAX2CoreNameSpaces, false );
    // AgentReader that reads genome and network
    NetAgentReader *doc_handler = new NetAgentReader( this, type );
    ErrorHandler *error_handler = (ErrorHandler*) doc_handler;
    parser->setContentHandler( doc_handler );
    parser->setErrorHandler( error_handler );    

    try {
        parser->parse( fname.c_str() );
    } catch( const XMLException &toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n" << message << "\n";
        XMLString::release( &message );
    } catch( const SAXParseException &toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n" << message << "\n";
        XMLString::release( &message );
    } catch( boost::bad_lexical_cast &toCatch ) {
        std::cout << "Bad lexical cast: " << toCatch.what() << "\n";
    } catch( ... ) {
        std::cout << "Unexpected Exception \n" ;
    }

    // get that agent
    GenRegAgent *result = 0;
    if( !doc_handler->sawErrors() ) {
        result = doc_handler->agent();
    }
    doc_handler->resetDocument();
    // delete parser before calling Terminate
    delete parser;
    delete doc_handler;
    // and call Terminate
    XMLPlatformUtils::Terminate();
    
    // Now we need to read the network from a dot file...
    return result;
}

boost::tuple< std::vector< blowhole::Agent* >, 
    std::vector< blowhole::location > >
blowhole::Factory::readNetPop( std::string fname, int tt ) {
    // Read genomes
    try {
        XMLPlatformUtils::Initialize();
    } catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Error during initialization! :\n" << message << "\n";
        XMLString::release( &message );
    }
            
    SAX2XMLReader *parser = XMLReaderFactory::createXMLReader();
    parser->setFeature( XMLUni::fgSAX2CoreValidation, false );
    parser->setFeature( XMLUni::fgSAX2CoreNameSpaces, false );

    // not giving type, coz it's read from file
    NetPopReader *doc_handler = new NetPopReader( this );
    ErrorHandler *error_handler = (ErrorHandler*) doc_handler;
    parser->setContentHandler( doc_handler );
    parser->setErrorHandler( error_handler );    

    try {
        parser->parse( fname.c_str() );
    } catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n" << message << "\n";
        XMLString::release( &message );
    } catch( const SAXParseException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n" << message << "\n";
        XMLString::release( &message );
    } catch( boost::bad_lexical_cast &toCatch ) {
        std::cout << "Bad lexical cast: " << toCatch.what() << "\n";
    } catch( ... ) {
        std::cout << "Unexpected Exception \n" ;
    }

    // get that agent population
    std::vector< Agent* > result;
    std::vector< location > loc;
    if( !doc_handler->sawErrors() ) {
        boost::tie( result, loc ) = doc_handler->population();
    }
    doc_handler->resetDocument();

    // delete parser before calling Terminate
    delete parser;
    delete doc_handler;
    // and call Terminate
    XMLPlatformUtils::Terminate();

    // still need to read networks from dot file
    network( tt );
    return boost::make_tuple( result, loc );
}

void
blowhole::Factory::readAgentConfigurations() {
    for( int i = 0; i != conf_->optionAsInt( "nr_agent_type" ); ++i ) {
        std::string fname =
            "agent" + boost::lexical_cast< std::string >( i ) + ".cfg";
        conf_->parseAgentFile( fname );
        setConfiguration( i );
    }
}

void
blowhole::Factory::setConfiguration( int type ) {
    // pre: conf_ is initialised
    SimpleAgent::birthRate( conf_->optionAsDouble( "birth_rate", type ) );
    SimpleAgent::deathRate( conf_->optionAsDouble( "death_rate", type ) );
    
    NetAgent::deathRate( conf_->optionAsDouble( "death_rate", type ) );
    NetAgent::maxDistance( conf_->optionAsInt( "max_distance", type ) );
    NetAgent::maxGenomeSize( conf_->optionAsInt( "max_genome_size", type ) );
    NetAgent::maxTposons( conf_->optionAsInt( "max_tposons", type ) );
    NetAgent::genomePenaltyRate( 
        conf_->optionAsDouble( "genome_size_penalty", type ) );
    NetAgent::retroposonPenaltyRate(
        conf_->optionAsDouble( "tposons_penalty", type ) );
    // read in initial state
    boost::dynamic_bitset<> bux( conf_->optionAsString( "init_state", type ) );
    NetAgent::initialState( bux );

    DeltaAgent::deathRate( conf_->optionAsDouble( "death_rate", type ) );
    DeltaAgent::maxDistance( conf_->optionAsInt( "max_distance", type ) );
    DeltaAgent::maxGenomeSize( conf_->optionAsInt( "max_genome_size", type ) );
    DeltaAgent::maxTposons( conf_->optionAsInt( "max_tposons", type ) );
    DeltaAgent::genomePenaltyRate( 
        conf_->optionAsDouble( "genome_size_penalty", type ) );
    DeltaAgent::retroposonPenaltyRate(
        conf_->optionAsDouble( "tposons_penalty", type ) );
    // read in initial state
    boost::dynamic_bitset<> cux( conf_->optionAsString( "init_state", type ) );
    DeltaAgent::initialState( cux );
    DeltaAgent::nrGenes( conf_->optionAsInt( "genes", type ) );
    
    Chromosome::copyGeneRate( conf_->optionAsDouble( "cp_gene", type ) );
    Chromosome::removeGeneRate( conf_->optionAsDouble( "rm_gene", type ) );
    Chromosome::copyRetroposonRate( conf_->optionAsDouble( "cp_tp", type ) );
    Chromosome::removeRetroposonRate( conf_->optionAsDouble( "rm_tp", type ) );
    Chromosome::newRetroposonRate( conf_->optionAsDouble( "new_tp", type ) );
    Chromosome::removeRepeatRate( conf_->optionAsDouble( "rm_ltr", type ) );
    Chromosome::recombinationRate( 
        conf_->optionAsDouble( "dsb_recombination", type ) );
    Chromosome::thresholdRate( conf_->optionAsDouble( "thr_rate", type ) );
    Chromosome::weightInteractionRate( 
        conf_->optionAsDouble( "weight_ia_rate", type ) );
    Chromosome::copyInteractionRate( 
        conf_->optionAsDouble( "cp_ia_rate", type ) );
    Chromosome::removeInteractionRate( 
        conf_->optionAsDouble( "rm_ia_rate", type ) );
    Chromosome::newInteractionRate( 
        conf_->optionAsDouble( "new_ia_rate", type ) );
    Chromosome::tagInteractionRate( 
        conf_->optionAsDouble( "tag_ia_rate", type ) );
}

std::vector< std::list< blowhole::ChromosomeElement* > > *
blowhole::Factory::chromosomeParts( int type ) {
    // # genes
    int n_ds = conf_->optionAsInt( "genes", type );
    // # retroposons ( inbetween other genes )
    int n_tp = conf_->optionAsInt( "tposons", type );
    // # single repeat elements ( inbetween genes etc )
    int n_rp = conf_->optionAsInt( "repeats", type );

    // create chromosome parts
    std::vector< std::list< ChromosomeElement* > > *result = 
        new std::vector< std::list< ChromosomeElement* > >();
    // create genes
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_ds; ++i ) {
        result->back().push_back( gene() );
    }
    // create retroposons
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_tp; ++i ) {
        result->back().push_back( retroposon() );
    }
    // create single repeats
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_rp; ++i ) {
        result->back().push_back( repeat() );
    }
    return result;
}

std::list< blowhole::ChromosomeElement* > *
blowhole::Factory::organiseChromosome( 
    std::vector< std::list< ChromosomeElement* > > *parts, double organise ) {
    typedef std::vector< std::list< ChromosomeElement* > >::iterator aux_iter;
    std::list< ChromosomeElement* > *result =
        new std::list< ChromosomeElement* >();
    
    // paste perfect & shuffle a bit
    // last 2 elements of 'parts' are retroposon and repeat elements
    std::list< ChromosomeElement* > ltr = parts->back();
    parts->pop_back();
    std::list< ChromosomeElement* > rp = parts->back();
    parts->pop_back();
    // add repeats
    int aux = 2 * parts->size();
    Chromosome::ce_iter i = ltr.begin();
    while( i != ltr.end() ) {
        int bux = static_cast< int >( aux * uniform() );
        if( bux % 2 == 0 ) {
            ( *parts )[ bux / 2 ].push_back( *i );
        } else {
            ( *parts )[ bux / 2 ].push_front( *i );
        }
        ++i;
    }
    // and retroposons
    i = rp.begin();
    while( i != rp.end() ) {
        int bux = static_cast< int >( aux * uniform() );
        if( bux % 2 == 0 ) {
            ( *parts )[ bux / 2 ].push_back( *i );
        } else {
            ( *parts )[ bux / 2 ].push_front( *i );
        }
        ++i;
    }
    // glue together
    for( aux_iter i = parts->begin(); i != parts->end(); ++i ) {
        result->splice( result->end(), *i );
    }
    // shuffle a little bit...
    aux = result->size();
    int nr_elem = static_cast< int >( 0.5 + ( 1.0 - organise ) * aux );
    if(  nr_elem == aux ) {
        // most annoying, copy to vector, from vector
        std::vector< ChromosomeElement* > bux( result->begin(), result->end() );
        std::random_shuffle( bux.begin(), bux.end(), rand_range< int > );
        std::copy( bux.begin(), bux.end(), result->begin() );
    } else {
        // naive implementation
        // take out elements
        std::list< ChromosomeElement* > bux;
        for( int i = 0; i < nr_elem; ++i ) {
            Chromosome::ce_iter cux = boost::next( result->begin(), 
                static_cast< int >( aux * uniform() ) );
            bux.push_back( *cux );
            result->erase( cux );
            --aux;
        }
        // put them back 
        while( !bux.empty() ) {
            result->insert( 
                boost::next( result->begin(), 
                    static_cast< int >( aux * uniform() ) ), bux.back() );
            bux.pop_back();
            ++aux;
        }
    }
    delete parts;
    return result;
}
