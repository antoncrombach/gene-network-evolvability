//
// XML sax reader for agents
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "defs.hh"
#include "factory.hh"

#include "net_agent_reader.hh"

#include "hop_graph.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "retroposon.hh"
#include "repeat.hh"

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif

blowhole::NetAgentReader::NetAgentReader( Factory *ft, int type ) 
    : factory_( ft ), 
      sawErrors_( false ), done_( false ), class_( false ), 
      ref_state_( false ), agent_class_(), type_( type ), 
      agent_( 0 ), network_( new HopfieldGraph() ), genome_( 0 ), chromo_( 0 ), 
      chr_( new std::list< ChromosomeElement* >() ) {}

void
blowhole::NetAgentReader::startElement( const XMLCh* const uri, 
        const XMLCh* const localname, const XMLCh* const qname,
        const Attributes &attrs ) {
    // now check type of element
    char *xxx = XMLString::transcode( localname );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( aux == "agent" ) {
        // create tag
        char *bux = XMLString::transcode( attrs.getValue( 0u ) );
        me_.time = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        me_.loc.first = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        me_.loc.second = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        
        bux = XMLString::transcode( attrs.getValue( 3 ) );
        me_.i = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
    } else if( aux == "class" ) {
        // we know what kind of agent to create
        class_ = true;
    } else if( aux == "refstate" ) {
        ref_state_ = true;
    } else if( aux == "parent" ) {
        // create parent tag
        char* bux = XMLString::transcode( attrs.getValue( 0u ) );
        mother_.time = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        mother_.loc.first = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        mother_.loc.second = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        
        bux = XMLString::transcode( attrs.getValue( 3 ) );
        mother_.i = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
    } else if( aux == "gene" ) {
        // create another gene
        // id
        char* bux = XMLString::transcode( attrs.getValue( 0u ) );
        int cux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // state
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        int dux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // threshold
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        int eux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // create gene
        chr_->push_back( new Gene( cux, dux, eux ) );
        // and count it
        gene_set_.insert( cux );
    } else if( aux == "repeat" ) {
        // create repeat
        chr_->push_back( new Repeat() );
    } else if( aux == "tposon" ) {
        // create tposon
        char* bux = XMLString::transcode( attrs.getValue( 0u ) );
        int cux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        chr_->push_back( new Retroposon( cux ) );
    } else if( aux == "ia" ) {
        // create interaction
        // target
        char* bux = XMLString::transcode( attrs.getValue( 0u ) );
        int cux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        // interaction value
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        int dux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );        
        chr_->push_back( new Interaction( cux, dux ) );
    }
}

void
blowhole::NetAgentReader::characters( const XMLCh* const chars, 
        const uint length) {
    char *xxx = XMLString::transcode( chars );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( class_ ) {
        agent_class_ = aux;
        class_ = false;
    }
    if( ref_state_ ) {
        reference_ = boost::dynamic_bitset<>( aux.erase( 0, 2 ) );
        ref_state_ = false;
    }
}

void
blowhole::NetAgentReader::endElement( const XMLCh* const uri, 
        const XMLCh* const localname, const XMLCh* const qname ) {
    char *xxx = XMLString::transcode( localname );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( aux == "agent" ) {
        // ready to create the right agent
        if( agent_class_ == "NetAgent" ) {
            agent_ = new NetAgent( type_, genome_ );
            agent_->myTag( me_ );
            agent_->parentTag( mother_ );
            agent_->readFromFile( true );
        } else if( agent_class_ == "DeltaAgent" ) {
            agent_ = new DeltaAgent( type_, genome_, reference_ );
            agent_->myTag( me_ );
            agent_->parentTag( mother_ );
            agent_->readFromFile( true );
        } else {
            throw "Unknown agent.";
        }
        done_ = true;
    } else if( aux == "genome" ) {
        // create genome (with only one chromosome)
        std::list< Chromosome* > *ll = new std::list< Chromosome* >();
        ll->push_back( chromo_ );
        genome_ = new Genome( ll );
    } else if( aux == "chromosome" ) {
        // create chromosome
        const Config &cf = factory_->configuration();
        chromo_ = new Chromosome( 0, chr_ );
        chromo_->copyGeneRate( cf.optionAsDouble( "cp_gene", type_ ) );
        chromo_->removeGeneRate( cf.optionAsDouble( "rm_gene", type_ ) );
        chromo_->copyRetroposonRate( cf.optionAsDouble( "cp_tp", type_ ) );
        chromo_->removeRetroposonRate( cf.optionAsDouble( "rm_tp", type_ ) );
        chromo_->newRetroposonRate( cf.optionAsDouble( "new_tp", type_ ) );
        chromo_->removeRepeatRate( cf.optionAsDouble( "rm_ltr", type_ ) );
        chromo_->recombinationRate( 
            cf.optionAsDouble( "dsb_recombination", type_ ) );
        chromo_->thresholdRate( cf.optionAsDouble( "thr_rate", type_ ) );
        chromo_->weightInteractionRate( 
            cf.optionAsDouble( "weight_ia_rate", type_ ) );
        chromo_->copyInteractionRate( 
            cf.optionAsDouble( "cp_ia_rate", type_ ) );
        chromo_->removeInteractionRate( 
            cf.optionAsDouble( "rm_ia_rate", type_ ) );
        chromo_->newInteractionRate( 
            cf.optionAsDouble( "new_ia_rate", type_ ) );
        chromo_->tagInteractionRate( 
            cf.optionAsDouble( "tag_ia_rate", type_ ) );
        HopfieldGraph::perturbRate( 
            cf.optionAsDouble( "perturb_rate", type_ ) );
        // extra stuff
        if( agent_class_ == "DeltaAgent" ) {
            DeltaAgent::nrGenes( gene_set_.size() );
        }
    }
}

void
blowhole::NetAgentReader::resetDocument() {
    // just set to zero, we do not own this agent
    chr_ = 0;
    chromo_ = 0;
    genome_ = 0;
    network_ = 0;
    agent_ = 0;
}

