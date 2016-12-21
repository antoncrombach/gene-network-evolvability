//
// Implementation of three environments.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "environment.hh"

blowhole::Environment::Environment() 
    : model_( 0 ), generator_env_( 18 ), uniform_env_( generator_env_ ),
    changed_( false ) {}

void
blowhole::Environment::initialise( uint seed ) {
    generator_env_.seed( seed );
    uniform_gen_type aux( generator_env_ );
    uniform_env_ = aux;
    
    if( model_ != 0 ) {
        model_->population().evaluate( *this );
    }
}

void
blowhole::Environment::finish() {
    notify();
}

//
// Constant environment
//
blowhole::ConstantEnvironment::ConstantEnvironment( int k ) 
    : Environment(), attractor_( 0 ) {}

//
// Periodic environment
//
blowhole::PeriodicEnvironment::PeriodicEnvironment( int k ) 
    : Environment(), attractor_index_( 0 ), period_( 0 ), offset_( 0 ),
      states_( k, 0 ) {}

void
blowhole::PeriodicEnvironment::initialise( uint seed ) {
    // ignoring s
    Environment::initialise( seed );
}

void
blowhole::PeriodicEnvironment::attractor( int k, 
    const boost::dynamic_bitset<> &low ) {
    states_[ k ] = low;
}

void
blowhole::PeriodicEnvironment::fluctuate( long time ) {
    bool aux = false;
    if( ( time - offset_ ) % period_ == 0 ) {
        aux = true;
        /*attractor_index_ = 1 - attractor_index_;*/
        ++attractor_index_;
        if( attractor_index_ == states_.size() ) attractor_index_ = 0;
        // inform observers
        notify();
    }
    if( aux ) {
        model_->population().evaluate( *this );
    }
}

//
// Poisson environment
//
blowhole::PoissonEnvironment::PoissonEnvironment( int k ) 
    : Environment(), attractor_index_( 0 ), lambda_( 0.0 ),
      states_( k, 0 ) {}

void
blowhole::PoissonEnvironment::attractor( int k, 
    const boost::dynamic_bitset<> &low ) {
    states_[ k ] = low;
}

void
blowhole::PoissonEnvironment::initialise( uint seed ) {
    Environment::initialise( seed );
}

void
blowhole::PoissonEnvironment::fluctuate( long time ) {
    bool aux = false;
    double uu = uniform_env_();
    if( uu < lambda_ ) {
        aux = true;
        /*attractor_index_ = 1 - attractor_index_;*/
        uint old_index_ = attractor_index_;
        do {
            attractor_index_ = 
                static_cast< int >( uniform_env_() * states_.size() );
        } while( attractor_index_ == old_index_ );
        // inform observers
        notify();
    }
    if( aux ) {
        changed_ = true;
        model_->population().evaluate( *this );
        changed_ = false;
    }
}

