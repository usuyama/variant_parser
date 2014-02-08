#!/usr/bin/env ruby

require 'pattern-match'

class Variant
  attr_accessor :type, :ref, :obs, :left, :right, :len

  # init by giving a hash
  def initialize(attributes={})
    attributes.each do |name, value|
      send("#{name}=", value)
    end
  end

  def self.from_cigar(c)
    v = Variant.new({:len => c[0..-2].to_i})
    case c[-1]
    when "M" then
      v.type = :match
    when "I" then
      v.type = :ins
      v.ref = "-"
    when "D" then
      v.type = :del
      v.obs = "-"
    when "S" then
      v.type = :soft
    else
      raise "not_supported_cigar"
    end
    v
  end

  def inspect
    "#{@type} #{@left}-#{@right} #{@ref} #{@obs} len:#{@len}"
  end
end

class MD
  attr_accessor :type, :len, :ref

  def initialize(tag)
    if tag[0] == "^"
      @type = :del
      @len = tag.length - 1
      @ref = tag[1..-1]
    elsif tag =~ /\d+/
      @type = :match
      @len = tag.to_i
    else
      @type = :mismatch
      @len = 1
      @ref = tag[0]
    end
  end
end

class Read
  attr_accessor :seq_name, :flag, :chr, :pos, :map_quality, :cigar, :seq, :base_quality, :md, :left_most_pos

  def initialize(r)
    @seq_name, @flag, @chr, @pos, @map_quality, @cigar = r[0..5]
    @pos = @pos.to_i
    @seq = r[9]
    @base_quality = r[10]
    @md = r[12]

    # pos does not include the soft-clipped bases
    cigar = @cigar.scan(/\d+./)
    if cigar[0][-1] == "S"
      @left_most_pos = @pos - cigar[0][0..-2].to_i
    else
      @left_most_pos = @pos
    end
  end

  def get_variants
    md = @md.scan(/\d+|[A-Z]+|\^[A-Z]+/)[2..-1].map {|m| MD.new(m) }

    # first, make variants from CIGAR
    cigar = @cigar.scan(/\d+./).map {|c| Variant.from_cigar(c) }

    # prepare index for tracking CIGAR & MD
    cigar_idx, md_idx = -1, -1
    current_cigar, current_md = nil, nil
    get_next_cigar = proc {
      cigar_idx += 1
      if cigar_idx < cigar.length
        current_cigar = cigar[cigar_idx]
      end
    }
    get_next_md = proc {
      md_idx += 1
      if md_idx < md.length
        current_md = md[md_idx]
      end
    }
    # prepare the first pair of CIGAR and MD
    get_next_cigar.call
    get_next_md.call

    @variants = []
    while cigar_idx < cigar.length and md_idx < md.length
      match([current_cigar.type, current_md.type]) do
        with _[:soft, _] {
          @variants << current_cigar.clone
          get_next_cigar.call
        }
        with _[_, :del] {
          current_cigar.ref = current_md.ref
          @variants << current_cigar.clone
          get_next_cigar.call
          get_next_md.call
        }
        with _[_, :mismatch] {
          @variants << Variant.new({:type => :mismatch, :ref => current_md.ref, :len => 1})
          current_cigar.len -= 1
          get_next_md.call
        }
        with _[:ins, :match] {
          @variants << current_cigar.clone
          get_next_cigar.call
        }
        with _[:match, :match] {
          if current_cigar.len > current_md.len
            @variants << Variant.new({:type => :match, :len => current_md.len })
            current_cigar.len -= current_md.len
            get_next_md.call
            get_next_cigar.call if current_cigar.len == 0
          else
            @variants << current_cigar.clone
            current_md.len -= current_cigar.len
            get_next_cigar.call
            get_next_md.call if current_md.len == 0
          end
        }
      end
    end

    # setting left & right pos for variants
    k = @left_most_pos
    for i in 0..@variants.length-1
      v = @variants[i]
      v.left = k - 1
      if v.type == :ins
        v.right = v.left
      else #v.type == :mismatch, :ins, :soft, :match
        k += v.len
        v.right = v.left + v.len
      end
    end

    # setting obs for variants
    b = 0
    for i in 0..@variants.length-1
      v = @variants[i]
      if v.type != :del
        v.obs = seq[b, v.len]
        b += v.len
      end
    end

    @variants
  end

end
