//
// Implementation of a few methods in the subject/observer pattern.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "observer.hh"
#include "stream_manager.hh"


blowhole::LogObserver::LogObserver( long i )
    : Observer(), interval_( i ), val_( 1 ) {
    log_ = 0;
}

blowhole::LogObserver::LogObserver( const LogObserver &lo ) 
    : Observer( lo ) {
    log_ = lo.log_;
    interval_ = lo.interval_;
    val_ = lo.val_;
}

blowhole::LogObserver::~LogObserver() {
    closeLog();
}

void 
blowhole::LogObserver::openLog( std::string fname ) {
    log_ = StreamManager::instance()->openOutFileStream( fname, 
            std::fstream::out | std::fstream::app );
}

void 
blowhole::LogObserver::closeLog() {
    if( log_ != 0 ) {
        finalize();
        StreamManager::instance()->closeOutFileStream( log_ );
        log_ = 0;
    }
}

void
blowhole::LogObserver::update( Subject *s ) {
    // template method
    if( val_ == 1 ) {
        val_ = interval_;
        doUpdate( s );
    } else {
        --val_;
    }
}

//
// async log observer
//
blowhole::AsyncLogObserver::AsyncLogObserver() : Observer() {
    log_ = 0;
}

blowhole::AsyncLogObserver::AsyncLogObserver( const AsyncLogObserver &lo ) 
    : Observer( lo ) {
    log_ = lo.log_;
}

blowhole::AsyncLogObserver::~AsyncLogObserver() {
    closeLog();
}

void 
blowhole::AsyncLogObserver::openLog( std::string fname ) {
    log_ = StreamManager::instance()->openOutFileStream( fname, 
            std::fstream::out | std::fstream::app );
}

void 
blowhole::AsyncLogObserver::closeLog() {
    if( log_ != 0 ) {
        finalize();
        StreamManager::instance()->closeOutFileStream( log_ );
        log_ = 0;
    }
}

