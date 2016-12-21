//
// Part of the observer/subject pattern
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "subject.hh"
#include "observer.hh"
#include "stream_manager.hh"

void 
blowhole::Subject::attach( Observer *o ) {
    obs_.push_back( o );
}

void 
blowhole::Subject::detach( Observer *o ) {
    // observer 'o' exists in obs_
    obs_.erase( std::find( obs_.begin(), obs_.end(), o ) );
}

void 
blowhole::Subject::detachAll() {
    obs_.clear();
}

void 
blowhole::Subject::notify() 
{ notify( this ); }

void
blowhole::Subject::notify( Subject *s ) {
    for( obs_iter i = obs_.begin(); i != obs_.end(); ++i ) {
        ( **i ).update( s );
    }
}

