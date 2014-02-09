#!/usr/bin/env ruby

require 'minitest/autorun'
require './sam_parser.rb'

$data = []
open "test.sam" do |f|
  while l = f.gets
    r = l.chomp.split("\t")
    $data.push Read.new(r)
  end
end

def should_include_variant(variants, type, ref, obs, left, right)
  assert((variants.find_index {|v|
    v.type == type and
    (type == :soft or
     (v.ref == ref and v.obs == obs)) and
    v.left == left and
    v.right == right
  }), [variants, type, ref, obs, left, right].inspect)
end

class TestRead < MiniTest::Unit::TestCase
  def test_basic
    r = $data[0]
    variants = r.get_variants
    should_include_variant(variants, :mismatch, "A", "T", 1824, 1825)
    should_include_variant(variants, :mismatch, "C", "A", 1846, 1847)
  end

  def test_del
    r = $data[1]
    variants = r.get_variants
    should_include_variant(variants, :del, "GAC", "-", 1826, 1829)
  end

  def test_ins
    r = $data[2]
    variants = r.get_variants
    should_include_variant(variants, :ins, "-", "GTC", 1328, 1328)
    should_include_variant(variants, :mismatch, "T", "G", 1338, 1339)
  end

  def test_softclip
    r = $data[3]
    variants = r.get_variants
    should_include_variant(variants, :soft, "GCAACAATGG", "T" * 10, 1368, 1378)
    should_include_variant(variants, :mismatch, "C", "A", 1419, 1420)
  end

  def test_ins2
    r = $data[4]
    variants = r.get_variants
    should_include_variant(variants, :mismatch, "C", "A", 1304, 1305)
    should_include_variant(variants, :ins, "-", "GTC", 1328, 1328)
  end

  def test_del2
    r = $data[5]
    variants = r.get_variants
    should_include_variant(variants, :mismatch, "T", "G", 1780, 1781)
    should_include_variant(variants, :del, "GAC", "-", 1826, 1829)
    should_include_variant(variants, :mismatch, "C", "A", 1846, 1847)
  end

end
