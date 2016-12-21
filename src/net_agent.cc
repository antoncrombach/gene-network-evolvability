//
// Implementation of a network agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "net_agent.hh"
#include "environment.hh"
#include "population.hh"


blowhole::NetAgent::NetAgent() : GenRegAgent() {}

blowhole::NetAgent::NetAgent( uint tag, Genome *g ) : GenRegAgent( tag, g ) {}

blowhole::NetAgent::NetAgent( const NetAgent &ag ) : GenRegAgent( ag ) {
    NetAgent::copy( ag );
}

blowhole::NetAgent::~NetAgent() {}

blowhole::Agent* 
blowhole::NetAgent::clone() const {
    return new NetAgent( *this );
}

void 
blowhole::NetAgent::copy( const Agent &ag ) {}

void 
blowhole::NetAgent::step( Population &pop ) {
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
        distance_ = network_->hamming( pop.model()->environment().attractor() );
    }
}

blowhole::Agent*
blowhole::NetAgent::sibling() {
    // change my parent info, coz after mitosis I'm not parent anymore
    parent_.size_ = genome_->length();
    parent_.distance_ = distance_;
    // first mitosis: now genome changes
    genome_->duplicate();
    genome_->mutate();
    Genome *sister_genome = genome_->split();
    // building sister net agent
    NetAgent *sister = new NetAgent( type_, sister_genome );
    // update parent info in sister
    sister->parent_.size_ = parent_.size_;
    sister->parent_.distance_ = parent_.distance_;
    return sister;
}

void
blowhole::NetAgent::evaluate( const Environment &env ) {
    bool lethal = false;
    if( network_->isEmpty() or genome_->hasMutation() ) {
        network_->build( *genome_ );
    }
    // successful build?
    if( not network_->hasMissingNode() ) {
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
        distance_ = network_->hamming( env.attractor() );
    }
}

void 
blowhole::NetAgent::writeXml( std::ostream &os ) const {
    // pre: genome_ != 0
    os << "<agent birth=\"" << me_.time 
       << "\" x=\"" << me_.loc.first << "\" y=\"" << me_.loc.second
       << "\" i=\"" << me_.i << "\">\n";
    os << "<parent birth=\"" << ancestor_.time
       << "\" x=\"" << ancestor_.loc.first << "\" y=\"" << ancestor_.loc.second
       << "\" i=\"" << ancestor_.i << "\"/>\n<class>NetAgent</class>\n";
    os << "<score>" << score() << "</score>\n<dist>" << distance_ 
       << "</dist>\n" << *genome_ << "</agent>\n";
}

