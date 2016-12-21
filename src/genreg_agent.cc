//
// Implementation of functions of any network agent.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "genreg_agent.hh"
#include "environment.hh"
#include "population.hh"


double blowhole::GenRegAgent::death_rate_ = 0.0;

int blowhole::GenRegAgent::max_dist_ = 0;
boost::dynamic_bitset<> blowhole::GenRegAgent::init_state_( 0 );
int blowhole::GenRegAgent::env_bit_ = -1;

int blowhole::GenRegAgent::penalty_genome_ = 0;
int blowhole::GenRegAgent::penalty_tposons_ = 0;
double blowhole::GenRegAgent::penalty_genome_rate_ = 0;
double blowhole::GenRegAgent::penalty_tposons_rate_ = 0;


blowhole::GenRegAgent::GenRegAgent() 
    : Agent( 0 ), score_( 2 * death_rate_ ), distance_( 0 ), 
      genome_( 0 ), network_( 0 ), read_from_file_( false ) {}

blowhole::GenRegAgent::GenRegAgent( uint tag ) 
    : Agent( tag ), score_( 2 * death_rate_ ), distance_( 0 ), 
      genome_( 0 ), network_( 0 ), read_from_file_( false ) {}
      
blowhole::GenRegAgent::GenRegAgent( uint tag, Genome *g ) 
    : Agent( tag ), score_( 0.0 ), distance_( 0 ), 
      genome_( g ), network_( new HopfieldGraph() ), 
      read_from_file_( false ) {}
    
blowhole::GenRegAgent::GenRegAgent( const GenRegAgent &ag ) 
    : Agent( 0 ), genome_( 0 ), network_( 0 ) {
    GenRegAgent::copy( ag );
}

blowhole::GenRegAgent::~GenRegAgent() {
    if( genome_ != 0 ) delete genome_;
    if( network_ != 0 ) delete network_;
}

void 
blowhole::GenRegAgent::copy( const Agent &ag ) {
    const GenRegAgent &sag = dynamic_cast< const GenRegAgent & >( ag );
    Agent::copy( ag );
    score_ = sag.score_;
    distance_ = sag.distance_;
    parent_.size_ = sag.parent_.size_;
    parent_.distance_ = sag.parent_.distance_;
    genome_ = sag.genome_->clone();
    // note: graph needs to be rebuilt!!
    network_ = sag.network_->clone();
    read_from_file_ = sag.read_from_file_;
}

void
blowhole::GenRegAgent::initialise() {
    // Make sure caching is ok
    // genome stuff
    genome_->nrRetroposons();
    genome_->length();
    if( !read_from_file_ ) {
        genome_->setState( init_state_ );
        genome_->initialise();
    }
    // network stuff.. not necessary 
}

double
blowhole::GenRegAgent::score() const { 
    double result = static_cast< double >( distance_ );
    // now see if we have any penalties to add
    int aux = genome_->length() - penalty_genome_;
    int bux = genome_->nrRetroposons() - penalty_tposons_;
    if( aux > 0 ) {
        result += aux * penalty_genome_rate_;
    }
    if( bux > 0 ) {
        result += bux * penalty_tposons_rate_;
    }
    // we've got the raw score now and calculate the fitness score
    if( result < max_dist_ ) {
        return 1.0 - ( result / static_cast< double >( max_dist_ ) );
    } else {
        return 0.0;
    }
}

void
blowhole::GenRegAgent::deathRate( double f )
{ death_rate_ = f; }

double
blowhole::GenRegAgent::deathRate()
{ return death_rate_; }

void
blowhole::GenRegAgent::maxDistance( int f ) 
{ max_dist_ = f; }

int
blowhole::GenRegAgent::maxDistance() 
{ return max_dist_; }

void
blowhole::GenRegAgent::maxGenomeSize( int f )
{ penalty_genome_ = f; }

void
blowhole::GenRegAgent::maxTposons( int f )
{ penalty_tposons_ = f; }

void
blowhole::GenRegAgent::genomePenaltyRate( double f )
{ penalty_genome_rate_ = f; }

void
blowhole::GenRegAgent::retroposonPenaltyRate( double f )
{ penalty_tposons_rate_ = f; }

void 
blowhole::GenRegAgent::initialState( const boost::dynamic_bitset<> &f )
{ init_state_ = f; }

void
blowhole::GenRegAgent::sensor( int f )
{ env_bit_ = f; }
