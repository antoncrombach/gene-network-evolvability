//
// Implementation of a special chromosome element, the repeat.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "repeat.hh"
#include "pool.hh"

template<> blowhole::ObjectCache< blowhole::Repeat >* 
blowhole::ObjectCache< blowhole::Repeat >::instance_ = 0;

blowhole::Repeat::Repeat() : ChromosomeElement(), dsb_( false ) {}

blowhole::Repeat::Repeat( const Repeat &tp ) : ChromosomeElement( tp ) {
    Repeat::copy( tp );
}

blowhole::ChromosomeElement* 
blowhole::Repeat::clone() const {
    //return new Repeat( *this );
    Repeat *ltr = ObjectCache< Repeat >::instance()->borrowObject();
    ltr->dsb_ = dsb_;
    return ltr;
}

void 
blowhole::Repeat::copy( const ChromosomeElement &ce ) {
    const Repeat *ltr = dynamic_cast< const Repeat * >( &ce );
    dsb_ = ltr->dsb_;
}

bool
blowhole::Repeat::toPool() {
    return ObjectCache< Repeat >::instance()->returnObject( this );
}

void
blowhole::Repeat::writeXml( std::ostream &os ) const { 
    os << "<repeat/>\n";
}

