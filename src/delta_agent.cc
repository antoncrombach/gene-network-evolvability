//
// Implementation of a network agent with a special fitness function.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "delta_agent.hh"
#include "environment.hh"
#include "population.hh"

int blowhole::DeltaAgent::nr_genes_ = 0;

blowhole::DeltaAgent::DeltaAgent() 
    : GenRegAgent(), ref_state_() {}

blowhole::DeltaAgent::DeltaAgent( uint tag, Genome *g, 
    const boost::dynamic_bitset<> &ref ) 
    : GenRegAgent( tag, g ), ref_state_( ref ) {}
    
blowhole::DeltaAgent::DeltaAgent( const DeltaAgent &ag ) : GenRegAgent( ag ) {
    DeltaAgent::copy( ag );
}

blowhole::DeltaAgent::~DeltaAgent() {}

blowhole::Agent* 
blowhole::DeltaAgent::clone() const {
    return new DeltaAgent( *this );
}

void 
blowhole::DeltaAgent::copy( const Agent &ag ) {
    const DeltaAgent &sag = dynamic_cast< const DeltaAgent & >( ag );
    ref_state_ = sag.ref_state_;
}

void 
blowhole::DeltaAgent::step( Population &pop ) {
    // are we dead?
    double rr = uniform();
    if( rr < death_rate_ or distance_ >= max_dist_ ) {
        dying_ = true;
    } 
    // we perturb a bit maybe and do a few propagation steps
    if( not dying_ ) {
        network_->perturb();
        // assuming network has been built
        if( not network_->inAttractor() ) {
            network_->propagate();
        }
        distance_ = nr_genes_ - network_->hamming( ref_state_ );
    }
}

blowhole::Agent*
blowhole::DeltaAgent::sibling() {
    // change my parent info, coz after mitosis I'm not parent anymore
    parent_.size_ = genome_->length();
    parent_.distance_ = distance_;
    
    // first mitosis
    genome_->duplicate();
    genome_->mutate();
    Genome *sister_genome = genome_->split();
    // building sister net agent
    DeltaAgent *sister = new DeltaAgent( type_, sister_genome, ref_state_ );
    
    // update parent info in sister
    sister->parent_.size_ = parent_.size_;
    sister->parent_.distance_ = parent_.distance_;
    return sister;
}

void
blowhole::DeltaAgent::evaluate( const Environment &env ) {
    bool lethal = false;
    if( network_->isEmpty() or genome_->hasMutation() ) {
        network_->build( *genome_ );
    }
    // successful build?
    if( !network_->hasMissingNode() ) {
#ifdef SENSOR
        // environmental sensor
        genome_->setState( env_bit_, env.attractor()[ env_bit_ ] );
        // sensor end
#endif
        network_->propagate();
        
    } else {
        lethal = true;
    }
    // what's the story? Check for all zeroes!
    if( lethal or network_->allZero() ) {
        distance_ = max_dist_;
    } else {
        if( env.hasChanged() ) updateReferenceState();
        distance_ = nr_genes_ - network_->hamming( ref_state_ );
    }
}

void
blowhole::DeltaAgent::updateReferenceState() {
    network_->state( ref_state_ );
}

void
blowhole::DeltaAgent::referenceState( const boost::dynamic_bitset<> &f ) {
    ref_state_ = f;
}    

void
blowhole::DeltaAgent::referenceState( std::vector< int > &f ) const {
    for( boost::dynamic_bitset<>::size_type i = 0; i < ref_state_.size(); ++i ){
        f.push_back( ref_state_[ i ] );
    }
}

void 
blowhole::DeltaAgent::writeXml( std::ostream &os ) const {
    // pre: genome_ != 0
    os << "<agent birth=\"" << me_.time 
       << "\" x=\"" << me_.loc.first << "\" y=\"" << me_.loc.second
       << "\" i=\"" << me_.i << "\">\n";
    os << "<parent birth=\"" << ancestor_.time
       << "\" x=\"" << ancestor_.loc.first << "\" y=\"" << ancestor_.loc.second
       << "\" i=\"" << ancestor_.i << "\"/>\n<class>DeltaAgent</class>\n";
    os << "<score>" << score() << "</score>\n<dist>" << distance_ 
       << "</dist>\n<refstate>0b" << ref_state_ << "</refstate>\n"
       << *genome_ << "</agent>\n";
}

void
blowhole::DeltaAgent::nrGenes( int f )
{ nr_genes_ = f; }

int
blowhole::DeltaAgent::nrGenes()
{ return nr_genes_; }
