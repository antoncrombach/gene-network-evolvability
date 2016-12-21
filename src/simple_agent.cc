//
// Implementation of a simple agent (not complete).
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "simple_agent.hh"
#include "population.hh"


double blowhole::SimpleAgent::birth_rate_ = 0.0;
double blowhole::SimpleAgent::death_rate_ = 0.0;


blowhole::SimpleAgent::SimpleAgent() : Agent( 0 ), score_( birth_rate_ ) {}

blowhole::SimpleAgent::SimpleAgent( const SimpleAgent &ag ) : Agent( 0 ) {
    SimpleAgent::copy( ag );
}

blowhole::Agent* 
blowhole::SimpleAgent::clone() const {
    return new SimpleAgent( *this );
}

void 
blowhole::SimpleAgent::copy( const Agent &ag ) {
    const SimpleAgent &sag = dynamic_cast< const SimpleAgent & >( ag );
    Agent::copy( ag );
    score_ = sag.score_;
}

void 
blowhole::SimpleAgent::step( Population &pop ) {
    double rr = uniform();
    if( rr < death_rate_ ) {
        dying_ = true;
    }
}

blowhole::Agent*
blowhole::SimpleAgent::sibling() {
    Agent* result;
    result = this->clone();
    //SimpleAgent *aux = dynamic_cast< SimpleAgent* >( result );
    return result;
}

void 
blowhole::SimpleAgent::writeXml( std::ostream &os ) const {
    os << "<agent birth=\"" << me_.time << "\" x=\"" << me_.loc.first;
    os << "\" y=\"" << me_.loc.second << "\" i=\"" << me_.i << "\"";
    os << ">\n<class>SimpleAgent</class>\n</agent>\n";
}

void
blowhole::SimpleAgent::birthRate( double f )
{ birth_rate_ = f; }

double
blowhole::SimpleAgent::birthRate()
{ return birth_rate_; }

void
blowhole::SimpleAgent::deathRate( double f )
{ death_rate_ = f; }

double
blowhole::SimpleAgent::deathRate()
{ return death_rate_; }

